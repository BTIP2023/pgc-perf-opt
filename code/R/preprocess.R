# File: preprocess.R
# Main function: preprocess
# Preprocessing is also done by the other steps (i.e. dim-reduce and clustering)
# so this is more accurately called `initial preprocess`.
# Using the raw GISAID data, this function performs:
# 1. Data extraction, parsing, and augmentation (i.e. wrangling)
# 2. Stratified random sampling
# 3. Data sanitation and grouping

# Auxiliary function(s):
# generate_interm: generate intermediate fasta and metadata files

# preprocess assumes .tar.gz filename format is "country-variant-...".
# Extract GISAID tars to data/GISAID/datasets/country-variant/
# If data/GISAID/datasets/ already exist, do not do this routine.
# data_path is GISAID data directory.
# extract_path is GISAID data extraction path after getting untarred.
# Note: Each tar = {tsv, fasta}
preprocess <- function(gisaid_data_path, gisaid_extract_path,
                       seed = 1234, strat_size = 100,
                       country_exposure = "Philippines",
                       write_fastacsv = FALSE, stamp) {
  # Extract GISAID data.
  if (dir.exists(gisaid_extract_path)) {
    message("GISAID data already extracted from tar archives.")
  } else {
    message("Extracting GISAID data to data/GISAID/datasets/...")
    tars <- list.files(gisaid_data_path, pattern = ".+\\.tar")
    for (file_name in tars) {
      subdir <- str_match(file_name, pattern = "[^-]+-[^-]+")
      message(paste("Extracting to:", subdir))
      untar(paste(gisaid_data_path, file_name, sep = "/"),
            exdir = paste(gisaid_extract_path, subdir, sep = "/"))
    }
  }
  
  # Merge all extracted fasta and tsv files.
  ## fasta contains sequence, while tsv contains metadata.
  omicron_sub = c("ba275", "xbb", "xbb_1.5", "xbb_1.16", "xbb1.91")
  
  fastas <- list.files(gisaid_extract_path, recursive = TRUE,
                       pattern = ".+\\.fasta")
  tsvs <- list.files(gisaid_extract_path, recursive = TRUE,
                     pattern = ".+\\.tsv")
  nfiles <- length(fastas)
  
  ## Accumulators: fasta_all and metadata_all
  fasta_all <- list()
  metadata_all <- tibble()

  ## Warnings suppressed for data parsing, but handled cleanly so don't worry.
  suppressWarnings({
    for (i in 1:nfiles) {
      fasta_path <- paste(gisaid_extract_path, fastas[i], sep = "/")
      tsv_path <- paste(gisaid_extract_path, tsvs[i], sep = "/")
      variant <- str_match(fasta_path, pattern = "(?<=-).*(?=\\/)")
      if (variant %vin% omicron_sub) {
        variant <- "Omicron Sub"
      }
      variant <- str_to_title(variant)
      
      message(paste0("\nReading ", fasta_path, "... "))
      # Parse then merge fasta file with accumulator.
      fasta <- ape::read.FASTA(fasta_path)
      fasta_all <- append(fasta_all, fasta)
      message("\bDONE.")
      
      message(paste0("Reading ", tsv_path, "... "))
      # Parse then merge metaData file with accumulator.
      # Defer sanitation after random sampling so fasta and metaData kept 1:1.
      # Can't directly col_types = "c_c_D____ccc___if_c__cc_____" because
      # of dirt in some columns, still need to use characters.
      metaData <- read_tsv(tsv_path,
                           col_select = c(1,3,5,10,11,12,16,17,19,22,23),
                           show_col_types = FALSE)
      message("\bDONE.")
      
      # Not removing raw date as I believe it is useful for sorting or can be
      # parsed on an as-needed basis. Dropped year, month, day: just extract
      # them from the date using lubridate::{year,month,day}(date).
      metaData <- metaData %>%
        dplyr::mutate(variant = as.character(variant))
      
      # Note: Cannot use tidyr::nest(fasta or tibble(fasta)), see reason below.
      
      # Coerce guessed column types to correct types (also for bind_rows).
      # NAs introduced by coercion will be dropped later because dropping
      # metadata rows must be consistent with dropping fasta entries for
      # the reason that nested DNAbin lists are not supported in R.
      metaData$age <- as.integer(metaData$age)
      metaData$sex <- as.character(metaData$sex)
      
      metadata_all <- bind_rows(metadata_all, metaData)
    }
  })

  rm(fasta)
  rm(metaData)

  # Print out total samples beforehand to guide future strat_size.
  message(paste("\nTotal number of samples in complete, unpruned data:",
              nrow(metadata_all)))

  # Addon: Filter by country_exposure.
  drop_idxs <- which(metadata_all$country_exposure != country_exposure)
  fasta_all <- fasta_all[is.na(pmatch(1:length(fasta_all), drop_idxs))]
  metadata_all <- metadata_all[is.na(pmatch(1:nrow(metadata_all), drop_idxs)),]
  
  rm(drop_idxs)
  
  # At this point, fasta_all and metadata_all contains the needed data.
  # Now do stratified random sampling.
  set.seed(seed)
  
  # Append rowname column for fasta sampling.
  # Drop this column before exporting.
  meta_grouped <- metadata_all %>%
    dplyr::group_by(variant) %>%
    tibble::rownames_to_column()
  
  # Do not preserve grouping structure (below only) to avoid NULL groups
  # If number of samples in variant group < strat_size, then filter from
  # meta_grouped and put temporarily in dropped_variants.
  # If number of samples is variant group >= strat_size, then get those
  # and place in meta_grouped, then randomly sample each of those groups
  dropped_variants <- filter(meta_grouped, n() < strat_size)
  meta_grouped <- filter(meta_grouped, n() >= strat_size)
  if (nrow(meta_grouped) >= strat_size)
    meta_grouped <- sample_n(meta_grouped, strat_size)
  metadata_all <- bind_rows(meta_grouped, dropped_variants)
  
  # Remove grouping information from tibble, let downstream handle it
  metadata_all <- dplyr::ungroup(metadata_all)
  
  rm(dropped_variants)
  rm(meta_grouped)
  
  set.seed(NULL)
  
  idxs <- as.integer(metadata_all$rowname)
  fasta_all <- fasta_all[idxs]
  
  # Drop rowname column
  metadata_all = subset(metadata_all, select = -c(rowname))  
  
  # Drop rows with NA values and type mismatches.
  # Get the idxs of the dropped metadata_all rows then drop them in fasta_all.
  drop_idxs1 <- which(is.na(metadata_all), arr.ind=TRUE)[,1]
  drop_idxs2 <- c(which(is.numeric(metadata_all$sex)),
                  which(!(metadata_all$sex %vin% list("Male", "Female"))))
  drop_idxs3 <- which(lengths(fasta_all) == 0)
  drop_idxs <- unique(c(drop_idxs1, drop_idxs2, drop_idxs3))
  
  # Dropping below is analogous to select inverse.
  # pmatch creates matches, val for match and NA for no match.
  # We only take those without matches, i.e. those that won't be dropped.
  fasta_all <- fasta_all[is.na(pmatch(1:length(fasta_all), drop_idxs))]
  metadata_all <- metadata_all[is.na(pmatch(1:nrow(metadata_all), drop_idxs)),]
  
  rm(idxs, drop_idxs1, drop_idxs2, drop_idxs3, drop_idxs)
  
  # At this point, data has been stratified and randomly sampled.
  
  message(paste("Number of randomly selected samples in stratified data:",
                nrow(metadata_all)))
  
  # Addon: Fix regions
  metadata_all$division_exposure <- case_match(
    metadata_all$division_exposure,
    "Bicol" ~ "Bicol Region",
    "Calabarzon" ~ "CALABARZON",
    "Mimaropa" ~ "MIMAROPA",
    "National Capital Region" ~ "NCR",
    "Cordillera Administrative Region" ~ "CAR",
    "Ilocos" ~ "Ilocos Region",
    "Davao" ~ "Davao Region",
    "Bangsamoro Autonomous Region in Muslim Mindanao" ~ "BARMM",
    "Autonomous Region In Muslim Mindanao(ARMM)" ~ "BARMM",
    "Autonomous Region In Muslim Mindanao" ~ "BARMM",
    "Soccsksargen" ~ "SOCCSKSARGEN",
    "Region III" ~ "Central Luzon",
    "Region IV-B" ~ "MIMAROPA",
    "Zamboanga" ~ "Zamboanga Peninsula",
    "Region IV-A" ~ "CALABARZON",
    "Region VIII (Eastern Visayas)" ~ "Eastern Visayas",
    "Region VIII" ~ "Eastern Visayas",
    "Region X (Northern Mindanao)" ~ "Northern Mindanao",
    "Region XI (Davao Region)" ~ "Davao Region",
    "Region XII (Soccsksargen)" ~ "SOCCSKSARGEN",
    "Metropolitan Manila" ~ "NCR",
    .default = metadata_all$division_exposure
  )
  
  # Addon: Add shortened Regions in Roman Numerals
  metadata_all <- metadata_all %>%
    dplyr::mutate(division_code = case_match(
      division_exposure,
      "Ilocos Region" ~ "I",
      "Cagayan Valley" ~ "II",
      "Central Luzon" ~ "III",
      "CALABARZON" ~ "IV-A",
      "MIMAROPA" ~ "IV-B",
      "Bicol Region" ~ "V",
      "Western Visayas" ~ "VI",
      "Central Visayas" ~ "VII",
      "Eastern Visayas" ~ "VIII",
      "Zamboanga Peninsula" ~ "IX",
      "Northern Mindanao" ~ "X",
      "Davao Region" ~ "XI",
      "SOCCSKSARGEN" ~ "XII",
      "Caraga" ~ "XIII",
      .default = division_exposure
    ), .after = division_exposure)
  
  # Addon: Add age_group, adjacent to age column
  # Age group reference: https://www.statcan.gc.ca/en/concepts/definitions/age2
  metadata_all <- metadata_all %>%
    dplyr::mutate(age_group = cut(age, breaks=c(0,14,24,64,500),
                                  include.lowest=T,
                                  labels=c("0-14", "15-24", "25-64", "65+")),
                  .after = age)
  
  # Return fasta_all and metadata_all, do final preprocess measures
  metadata_all$sex <- as.factor(metadata_all$sex)
  metadata_all$variant <- as.factor(metadata_all$variant)
  list(fasta_all, metadata_all)
}

generate_interm <- function(metadata_all, fasta_all) {
  message("\nWriting generated fasta and csv files:")
  message(paste0("Writing intermediate fasta to ",
                 sprintf("data/interm/fasta_all_%s.fasta... ", stamp)),
          appendLF = FALSE)
  ape::as.character.DNAbin(fasta_all)
  seqinr::write.fasta(fasta_all, names(fasta_all),
                      sprintf("data/interm/fasta_all_%s.fasta", stamp))
  message("DONE.")
  
  message(paste0("Writing intermediate metadata to ",
                 sprintf("data/interm/metadata_all_%s.csv... ", stamp)),
          appendLF = FALSE)
  write_csv(metadata_all,
            sprintf("data/interm/metadata_all_%s.csv", stamp))
  message("DONE.")
  
  # Now, get overview for the data that has been sampled.
  # Only sampled rows will be summarized.
  # compile_overview(metadata_all, 'data/overview')
  
  # After getting credits, we can now drop submitting_lab and authors
  metadata_all <- subset(metadata_all, select = -c(submitting_lab, authors))
  
  # Refetch fasta_all data using ape::read.FASTA to optimize for kcount.
  # kcount using DNAbin is faster than characters.
  fasta_all <- read.FASTA(sprintf("data/interm/fasta_all_%s.fasta", stamp))
}
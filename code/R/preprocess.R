# File: preprocess.R
# Main function: preprocess
# Note: See defaults in definition.
# Using the raw GISAID data, this function performs:
# 1. Data extraction, parsing, and augmentation
# 2. Stratified random sampling
# 3. Data Sanitation and Grouping
# 4. Generation of intermediate fasta and metadata files

# Assumption: .tar.gz filename format is "country-variant-...".
# Extract GISAID tars to data/GISAID/datasets/country-variant/
# If data/GISAID/datasets/ already exist, do not do this routine.
# data_path is GISAID data directory.
# extract_path is GISAID data extraction path after getting untarred.
# Note: Each tar = {tsv, fasta}

preprocess <- function(data_path, extract_path,
                       seed = 1234, strat_size = 100,
                       country_exposure = "Philippines",
                       write_fastacsv = FALSE, stamp) {
  
  # Extract GISAID data.
  if (dir.exists(extract_path)) {
    message("GISAID data already extracted from tar archives.")
  } else {
    message("Extracting GISAID data to data/GISAID/datasets/...")
    tars <- list.files(data_path, pattern = ".+\\.tar")
    for (file_name in tars) {
      subdir <- str_match(file_name, pattern = "[^-]+-[^-]+")
      message(paste("Extracting to:", subdir))
      untar(paste(data_path, file_name, sep = "/"),
            exdir = paste(extract_path, subdir, sep = "/"))
    }
  }
  
  # Merge all extracted fasta and tsv files.
  ## fasta contains sequence, while tsv contains metadata.
  omicron_sub = c("ba275", "xbb", "xbb_1.5", "xbb_1.16", "xbb1.91")
  
  fastas <- list.files(extract_path, recursive = TRUE, pattern = ".+\\.fasta")
  tsvs <- list.files(extract_path, recursive = TRUE, pattern = ".+\\.tsv")
  nfiles <- length(fastas)
  
  ## Initialize accumulator data frames (faster than tibbles).
  fasta_all <- data.frame()
  metadata_all <- data.frame()
  
  ## Warnings suppressed for data parsing, but handled cleanly so don't worry.
  suppressWarnings({
    for (i in 1:nfiles) {
      fasta_path <- paste(extract_path, fastas[i], sep = "/")
      tsv_path <- paste(extract_path, tsvs[i], sep = "/")
      variant <- str_match(fasta_path, pattern = "(?<=-).*(?=\\/)")
      if (variant %vin% omicron_sub) {
        variant <- "Omicron Sub"
      }
      variant <- str_to_title(variant)
      
      message(paste0("\nReading ", fasta_path, "... "))
      # Parse then merge fasta file with accumulator.
      # Optimization: If write_fasta == TRUE, then use seqinr, else use ape.
      if (write_fastacsv) {
        fasta <- seqinr::read.fasta(fasta_path, forceDNAtolower = FALSE)
      } else {
        fasta <- ape::read.FASTA(fasta_path)
      }
      fasta_all <- c(fasta_all, fasta)
      message("\bDONE.")
      
      message(paste0("Reading ", tsv_path, "... "))
      # Parse then merge metaData file with accumulator.
      # Defer sanitation after random sampling so fasta and metaData kept 1:1.
      metaData <- as.data.frame(read_tsv(tsv_path,
                                         col_select = c(1,3,5,10,11,
                                                        12,16,17,19,22,23),
                                         show_col_types = FALSE))
      message("\bDONE.")
      
      # Not removing raw date as I believe it is useful for sorting or can be
      # parsed on an as-needed basis. Dropped year, month, day: just extract
      # them from the date using lubridate::{year,month,day}(date).
      metaData <- metaData %>%
        dplyr::mutate(variant = as.character(variant))
      
      # Coerce some columns to correct types.
      # NAs introduced by coercion will be dropped later because dropping
      # metadata rows must be consistent with dropping fasta entries for
      # the reason that efficient df columns with lists are not well-supported
      # in R. Tibbles on the other hand reduce efficiency and compatibility.
      # Also, using tibbles add another need to extract the fasta from the tibble
      # for later use (and for writeback), and rowwise operations in tibbles
      # are said to be slow.
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
  drop_idxs <- which(metadata_all$country != country_exposure)
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
  
  # At this point, data has been stratified and randomly sampled.
  # Now, get credits for the data that has been sampled.
  # Only sampled rows will be credited.
  # compile_overview(metadata_all)
  
  # After getting credits, we can now drop submitting_lab and authors
  metadata_all <- subset(metadata_all, select = -c(submitting_lab, authors))
  
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
  
  # Lines below creates intermediate fasta_all.fasta and metadata_all.csv.
  # Optimization: Check job order if want to write fasta and csv.
  if (write_fastacsv) {
    message("\nWriting generated fasta and csv files:")
    
    # Write parameters used to log file
    write_to_log(output_dir = "data/interm", filename = "log.txt",
                 log_string = sprintf("timestamp = %s\nseed = %d, strat_size = %d",
                                      stamp, seed, strat_size))
    
    message(paste0("Writing intermediate fasta to ",
                   sprintf("data/interm/fasta_all_%s.fasta... ", stamp)),
            appendLF = FALSE)
    
    seqinr::write.fasta(fasta_all, names(fasta_all),
                        sprintf("data/interm/fasta_all_%s.fasta", stamp))
    
    message("DONE.")
    
    message(paste0("Writing intermediate metadata to ",
                  sprintf("data/interm/metadata_all_%s.csv... ", stamp)),
            appendLF = FALSE)
    
    write.csv(metadata_all,
              sprintf("data/interm/metadata_all_%s.csv", stamp),
              row.names = FALSE)
    
    message("DONE.")
    
    # Refetch fasta_all data using ape::read.FASTA to optimize for kcount.
    fasta_all <- read.FASTA(sprintf("data/interm/fasta_all_%s.fasta", stamp))
  }
  
  # Return fasta_all and metadata_all
  list(fasta_all, metadata_all)
}
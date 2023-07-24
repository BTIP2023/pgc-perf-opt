# File: preprocess.R
# Preprocessing is also done by the other steps (i.e. dim-reduce and clustering)
# so this is more accurately called `initial preprocess`.

# Main function: get_sample
# Using the raw GISAID data, this function performs:
# 1. Data extraction, parsing, and augmentation (i.e. wrangling)
# 2. Stratified random sampling

# Auxiliary function(s):
# sanitize_sample: fix metadata column types and values
# generate_interm: generate intermediate fasta and metadata files
# compile_overview: generate overview of metadata, mainly summaries and credits
# generate_treemaps: generate treemaps to improve visualization of sampled data

# Note: All dropping of rows only done in get_sample.

# preprocess assumes .tar.gz filename format is "country-variant-...".
# Extract GISAID tars to data/GISAID/datasets/country-variant/
# If data/GISAID/datasets/ already exist, do not do this routine.
# data_path is GISAID data directory.
# extract_path is GISAID data extraction path after getting untarred.
# Note: Each tar = {tsv, fasta}
get_sample <- function(gisaid_data_path = "data/GISAID",
                       gisaid_extract_path = "data/GISAID/datasets",
                       seed = 1234, strat_size = 100,
                       country_exposure = "Philippines", stamp) {
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
      metaData <- readr::read_tsv(tsv_path,
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

  # Print out number of samples beforehand to guide future strat_size.
  message(paste("\nTotal number of samples in complete, unpruned data:",
              nrow(metadata_all)))
  message("Variant distribution in complete data:")
  metadata_all %>%
    dplyr::group_by(variant) %>%
    dplyr::count() %>% print()

  # Addon: Filter by country_exposure.
  drop_idxs <- which(metadata_all$country_exposure != country_exposure)
  fasta_all <- fasta_all[is.na(pmatch(1:length(fasta_all), drop_idxs))]
  metadata_all <- metadata_all[is.na(pmatch(1:nrow(metadata_all), drop_idxs)),]
  
  rm(drop_idxs)
  
  # At this point, fasta_all and metadata_all contains the needed data.
  # Now do stratified random sampling.
  set.seed(seed)
  
  # Append rowname column for fasta_all subsetting.
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
  # TODO: consider sample_frac()
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
  
  # Drop explicit rowname column, already used to subset fasta_all.
  metadata_all <- metadata_all %>% select(!rowname)
  
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
  
  # Addon: Add age_group, adjacent to age column
  metadata_all <- metadata_all %>%
    dplyr::mutate(age_group = cut(age, breaks=c(0,14,24,64,500),
                                  include.lowest=T,
                                  labels=c("0-14", "15-24", "25-64", "65+")),
                  .after = age)
  
  # At this point, data has been stratified and randomly sampled.
  # We may now return it for further cleaning and downstream use.
  message(paste("\nNumber of randomly selected samples in stratified data:",
                nrow(metadata_all)))
  message("Variant distribution in selected samples:")
  metadata_all %>%
    dplyr::group_by(variant) %>%
    dplyr::count() %>% print()
  
  list(fasta_all, metadata_all)
}

# Note that this assumes that country_exposure was set to "Philippines"
sanitize_sample <- function(metadata_all) {
  # Fix regions, following PH Atlas convention aside from MIMAROPA,
  # https://www.philatlas.com/regions.html
  message("Cleaning division_exposure and adding division_code... ",
          appendLF = FALSE)
  metadata_all <- metadata_all %>%
    dplyr::mutate(division_exposure = case_match(
      division_exposure,
      "Bicol" ~ "Bicol Region",
      "Calabarzon" ~ "CALABARZON",
      "Mimaropa" ~ "MIMAROPA",
      "Ilocos" ~ "Ilocos Region",
      "Davao" ~ "Davao Region",
      "Autonomous Region In Muslim Mindanao(ARMM)" ~
        "Bangsamoro Autonomous Region in Muslim Mindanao",
      "Autonomous Region In Muslim Mindanao" ~
        "Bangsamoro Autonomous Region in Muslim Mindanao",
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
      "Metropolitan Manila" ~ "National Capital Region",
      .default = metadata_all$division_exposure))
  
  # Add shortened Regions in Roman Numerals, call it division_code
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
      "Cordillera Administrative Region" ~ "CAR",
      "Bangsamoro Autonomous Region in Muslim Mindanao" ~ "BARMM",
      "National Capital Region" ~ "NCR",
      .default = division_exposure
    ), .after = division_exposure)
  message("DONE.")
  
  # Clean up submitting labs
  # Optimized .*(?=.* -) to .*(?= -.*)
  # Labs with "name1 - name2" will only retain name1
  # Labs with "/" not separated into rows, converted "/" to "and"
  message("Cleaning submitting_lab... ", appendLF = FALSE)
  metadata_all <- metadata_all %>%
    dplyr::mutate(submitting_lab = dplyr::case_when(
      stringr::str_detect(submitting_lab, regex("\\w/\\w")) ~
        stringr::str_replace(submitting_lab, "/", " and "),
      stringr::str_detect(submitting_lab, regex("Center Visayas",
                                                ignore_case = TRUE)) ~
        "Philippine Genome Center Visayas",
      stringr::str_detect(submitting_lab, "-") ~
        stringr::str_extract(submitting_lab, ".*(?= -.*)"),
      .default = submitting_lab
    ))
  message("DONE.")
  
  # Clean up authors
  # ",(?![A-Z]+)"
  message("Cleaning authors... ", appendLF = FALSE)
  metadata_all <- metadata_all %>%
    tidyr::separate_rows(authors, sep = ",| and |nE|. Chel") %>%
    dplyr::mutate(authors =
                    str_replace(authors, "Dr.|PhD|MSc|MD|RMT|FPSP", "")) %>%
    dplyr::mutate(authors = stringr::str_replace(authors, "√±", "ñ")) %>%
    dplyr::mutate(authors = stringr::str_squish(authors)) %>%
    dplyr::filter(authors != "") %>%
    dplyr::mutate(authors = dplyr::case_match(
      authors,
      "Eva Maria C. Cutiongco-de la" ~ "Eva Maria C. Cutiongco-de la Paz",
      "Eva Maria Cutiongco-de la Paz" ~ "Eva Maria C. Cutiongco-de la Paz",
      "Jefferson Earl Halog" ~ "Jefferson Earl J. Halog",
      "Joana Ina Manalo" ~ "Joanna Ina G. Manalo",
      "lcid Aaron R. Pangilinan" ~ "Elcid Aaron R. Pangilinan",
      "Marissa" ~ "Marissa M. Alejandria",
      "Nina Francesca Bustamante" ~ "Niña Francesca M. Bustamante",
      "Niña Francesca Bustamante" ~ "Niña Francesca M. Bustamante",
      "Renato Jacinto Q. Manta" ~ "Renato Jacinto Q. Mantaring",
      "Samantha Louise Bado" ~ "Samantha Luoise P. Bado",
      "Vanessa Joy Diamante" ~ "Vanessa Joy F. Diamante",
      "Yvonne Valerie Austria" ~ "Yvonne Valerie D. Austria",
      "sea Joy M. Galutan" ~ "Chelsea Joy M. Galutan",
      "Chelsea Joy Galutan" ~ "Chelsea Joy M. Galutan",
      "Celia Carlos" ~ "Celia C. Carlos",
      "Alyssa Joyce Telles" ~ "Alyssa Joyce E. Telles",
      "Ardiane Ysabelle Dolor" ~ "Ardiane Ysabelle M. Dolor",
      "Charalyn Babida" ~ "Charalyn A. Babida",
      "Cynthia P. S" ~ "Cynthia P. Saloma",
      "Devon Ray Pacial" ~ "Devon Ray O. Pacial",
      "Diadem Ricarte" ~ "Diadem R. Ricarte",
      "Edsel Maurice Salvana" ~ "Edsel Maurice T. Salvaña",
      "Edsel Maurice T. Salva" ~ "Edsel Maurice T. Salvaña",
      "Edsel Maurice T. Salvana" ~ "Edsel Maurice T. Salvaña",
      "Florrianne D. Collantes" ~ "Florianne D. Collantes",
      "Francisco Gerardo Polotan" ~ "Francisco Gerardo M. Polotan",
      "Gerald Ivan Sotelo" ~ "Gerald Ivan S. Sotelo",
      "Henrietta Marie Rodriguez" ~ "Henrietta Marie M. Rodriguez",
      "J oshua Gregor E. Dizon" ~ "Joshua Gregor E. Dizon",
      "Johnny A. Ong." ~ "Johnny A. Ong",
      "Alethea R. De Guzma" ~ "Alethea R. de Guzman",
      "Alethea R. de Guzma" ~ "Alethea R. de Guzman",
      "Alethea R. De Guzman" ~ "Alethea R. de Guzman",
      "Althea R. De Guzman" ~ "Alethea R. de Guzman",
      "Diomedes A. Carino" ~ "Diomedes A. Cariño",
      "Edsel Maurice Salvaña" ~ "Edsel Maurice T. Salvaña",
      "Anna Ong-Lim" ~ "Anna Lisa T. Ong-Lim",
      "Catalino Demetria" ~ "Catalino S. Demetria",
      "Fe C.Villarama" ~ "Fe C. Villarama",
      "John Michael Egana" ~ "John Michael C. Egana",
      "Joshua Jose Endozo"~ "Joshua Jose S. Endozo",
      "June Jay B.Tejano" ~ "June Jay B. Tejano",
      "Krisitna Patriz Dela Cruz" ~ "Kristina Patriz Dela Cruz",
      "Lei Lanna Dancel" ~ "Lei Lanna M. Dancel",
      "Lindsay Clare D.L. Carandang" ~ "Lindsay Claire D.L. Carandang",
      "Liza Mae De La Cruz" ~ "Liza Mae L. De La Cruz",
      "Ma Angelica Tujan" ~ "Ma. Angelica A. Tujan",
      "Ma Angelica A. Tujan" ~ "Ma. Angelica A. Tujan",
      "Ma. Exanil Plantig" ~ "Ma. Exanil L. Plantig",
      "Marielle M Gamboa" ~ "Marielle M. Gamboa",
      "Marissa Alejandria" ~ "Marissa M. Alejandria",
      "Renalyn SeNoran" ~ "Renalyn Señoran",
      "Rona Clarisse B. Sobreca" ~ "Rona Clarisse B. Sobrecarey",
      "Sophia Isabelle V. Diño" ~ "Sofia Isabelle V. Diño",
      "Anne M. eco" ~ "Anne M. Eco",
      "Elizabeth Freda O.Telan" ~ "Elizabeth Freda O. Telan",
      "Ma. Angelica Tujan" ~ "Ma. Angelica A. Tujan",
      "Mariko Siato-Obata" ~ "Mariko Saito-Obata",
      "April Mae Numeron" ~ "April Mae M. Numeron",
      .default = authors
    ))

  # Collapse authors back to authors column
  # TODO: This seems redundant as we will unravel this again
  # for compile_overview, so think of workaround.
  # Bug with .by in mutate, so use group_by before mutate.
  metadata_all <- metadata_all %>%
    dplyr::group_by(strain) %>%
    dplyr::arrange(authors) %>%
    dplyr::mutate(authors = paste(authors, collapse = ", ")) %>%
    dplyr::distinct(strain, variant, .keep_all = TRUE) %>%
    dplyr::ungroup()
  message("DONE.")
  
  # Add ph_region which is of the form "division_exposure (division_code)"
  message("Adding ph_region to metadata... ", appendLF = FALSE)
  abbregions <- list("BARMM", "CAR", "NCR")
  metadata_all <- metadata_all %>%
    dplyr::mutate(ph_region = dplyr::if_else(division_code %vin% abbregions,
      stringr::str_glue("{division_exposure} ({division_code})"),
      stringr::str_glue("{division_exposure} (Region {division_code})")),
      .after = division_code)
  message("DONE.")
  
  # Coerce sex and variant to factors to save space
  metadata_all$sex <- as.factor(metadata_all$sex)
  metadata_all$variant <- as.factor(metadata_all$variant)
  
  return(metadata_all)
}

# Now always writes intermediate files. Thanks ape.
generate_interm <- function(fasta_all, metadata_all,
                            write_path = "data/interm", stamp) {
  if (!dir.exists(write_path)) {
    dir.create(write_path)
  }
  fasta_path <- sprintf("%s/fasta_all_%s.fasta", write_path, stamp)
  csv_path <- sprintf("%s/metadata_all_%s.csv", write_path, stamp)
  
  message("Writing generated fasta and csv files:")
  message(sprintf("Writing intermediate fasta to %s... ", fasta_path),
          appendLF = FALSE)
  ape::write.FASTA(fasta_all, file = fasta_path)
  message("DONE.")
  message(sprintf("Writing intermediate metadata to %s... ", csv_path),
          appendLF = FALSE)
  readr::write_csv(metadata_all, file = csv_path)
  message(sprintf("Writing intermediate metadata to %s... DONE.", csv_path))
}

# Compile overview of sampled data
# Ex. Group by age_group and variant then count()
# Ex. Group by authors and how many samples they've submitted
# Only sampled rows will be given overviews.
compile_overview <- function(metadata_all, write_path = "data/overview") {
  # Get accession numbers and compile to a list
  gisaid_esp_isl <- sort(metadata_all$gisaid_epi_isl)
  
  # Get authors and number of samples they've submitted: per variant and total n
  df_authors <- metadata_all %>%
    tidyr::separate_rows(authors, sep = ", ") %>%
    dplyr::group_by(authors) %>%
    variants_per_factor() %>%
    dplyr::mutate(n = rowSums(across(where(is.numeric))))
  
  # Get submitting labs and the number of samples they have submitted:
  ## per variant and total n
  df_labs <- metadata_all %>%
    dplyr::group_by(submitting_lab) %>%
    variants_per_factor() %>%
    dplyr::mutate(n = rowSums(across(where(is.numeric))))
  
  # Get number of variants and total samples n per division_exposure
  # Included division_code and ph_region for utilitarian purposes
  df_division <- metadata_all %>%
    dplyr::group_by(division_exposure, division_code, ph_region) %>%
    variants_per_factor() %>%
    dplyr::mutate(n = rowSums(across(where(is.numeric))))
  
  # Get number of variants and total samples n per sex
  df_sex <- metadata_all %>%
    dplyr::group_by(sex) %>%
    variants_per_factor() %>%
    dplyr::mutate(n = rowSums(across(where(is.numeric))))
  
  # Get number of variants and total samples n per age_group
  df_age_group <- metadata_all %>%
    dplyr::group_by(age_group) %>%
    variants_per_factor() %>%
    dplyr::mutate(n = rowSums(across(where(is.numeric))))
  
  # WRITE overviews to write_path
  if (!dir.exists(write_path)) {
    dir.create(write_path)
  }
  
  # Write GISAID Accession Numbers
  write_lines(gisaid_esp_isl, paste(write_path, "accession.txt", sep = "/"))
  
  # Write the rest of overviews to CSVs
  write_csv(df_authors, paste(write_path, "authors.txt", sep = "/"))
  write_csv(df_labs, paste(write_path, "labs.csv", sep = "/"))
  write_csv(df_division, paste(write_path, "division.csv", sep = "/"))
  write_csv(df_age_group, paste(write_path, "age_group.csv", sep = "/"))
  write_csv(df_sex, paste(write_path, "sex.csv", sep = "/"))
  
  # After getting credits, we can now drop submitting_lab and authors
  metadata_all <- subset(metadata_all, select = -c(submitting_lab, authors))
  
  return(metadata_all)
}

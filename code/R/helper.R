# File: helper.R
# To contain auxiliary functions for all code/R/ files.

# Appends system information and parameters used to a text file.
# specified by output_path and filename.
# This function only supports Windows 10/11 and Linux.
# In Windows, assumes that there's only one processor installed.
write_to_log <- function(output_dir, filename, log_string, stamp) {
  if (!dir.exists(output_dir))
    dir.create(output_dir)
  output_path <- paste(output_dir, filename, sep = "/")
  
  fileConn<-file(output_path)
  
  if (pacman::p_detectOS() == "Windows") {
    systeminfo <- system("systeminfo", intern = TRUE)
    systeminfo <- systeminfo[c(3,4,13,14,15,16,17,25)]
    cpuinfo <- system("WMIC CPU Get DeviceID,NumberOfCores,NumberOfLogicalProcessors",
                      intern = TRUE)
    cpuinfo <- cpuinfo[-length(cpuinfo)]
    specs <- c(systeminfo, "----------------------------", cpuinfo)
    readr::write_lines(c(sprintf("=============%s=============", stamp),
                         specs, log_string),
                       file = output_path, append = TRUE)
  } else if (pacman::p_detectOS() == "Linux") {
    device <- paste(as.list(Sys.info())[c("sysname", "release")], collapse = ' ')
    lsb_release <- as.list(system("lsb_release -a", intern = TRUE,
                                  ignore.stderr = TRUE))
    names(lsb_release) <- str_match(lsb_release, pattern = ".+?(?=:)")
    lsb_release <- unlist(lsb_release[c("Description", "Codename")],
                          use.names = FALSE)
    lscpu <- as.list(system("lscpu", intern = TRUE))
    names(lscpu) <- str_match(lscpu, pattern = ".+?(?=:)")
    lscpu <- unlist(lscpu[c("Architecture", "CPU(s)", "Model name",
                            "Thread(s) per core", "Core(s) per socket",
                            "Virtualization", "Vulnerability Itlb multihit",
                            "Vulnerability L1tf",
                            "Vulnerability Mds",
                            "Vulnerability Meltdown",
                            "Vulnerability Mmio stale data",
                            "Vulnerability Retbleed",
                            "Vulnerability Spec store bypass",
                            "Vulnerability Spectre v1",
                            "Vulnerability Spectre v2",
                            "Vulnerability Srbds",
                            "Vulnerability Tsx async abort")],
                    use.names = FALSE)
    mem <- system("grep MemTotal /proc/meminfo", intern=T)
    specs <- c(device, lsb_release, lscpu, mem)
    readr::write_lines(c(sprintf("=============%s=============", stamp),
                         specs, log_string),
                       file = output_path, append = TRUE)
  } else {
    readr::write_lines(c(sprintf("=============%s=============", stamp),
                         "OS not supported by logger!", log_string),
                       file = output_path, append = TRUE)
  }
  close(fileConn)
}

# Returned time is of the format YYYYMMDDHHMMSSNNN
get_time <- function() {
  ret <- gsub("[.:-]|\\s", "", as.character(Sys.time()))
  ret <- substring(ret, 1, nchar(ret)-3)
}

# Wrapper function for making treemaps with variable ...length() levels
# NOTE: The treemap function can plot any treemap you can think of, yeah!
treemap <- function(df, ..., tm_title = "",
                    tm_subtitle = "", tm_caption = "") {
  # Create df containing the columns to summarize.
  summ <- df %>% select(...) %>% tibble(n = rep(1, nrow(df)))
  # Generate level JSONs
  lvl_opts <- list()
  for (i in 1:...length()) {
    if (i == 1) {
      lvl_opts[[i]] <- list(
        level = 1,
        borderWidth = 2,
        borderColor = "white",
        colorVariation = list(
          key = "brightness",
          to = 0.250
        ),
        dataLabels = list(
          enabled = TRUE,
          align = "left",
          verticalAlign = "top",
          style = list(
            fontSize = "1.2em",
            color = "white",
            textOutline = FALSE,
            fontWeight = "normal")))
    } else if (i == 2) {
      lvl_opts[[i]] <- list(
        level = 2,
        colorVariation = list(
          key = "brightness",
          to = 0.250
        ),
        dataLabels = list(
          enabled = TRUE,
          style = list(
            fontSize = "0.8em",
            color = "white",
            textOutline = FALSE,
            fontWeight = "normal")))
    } else {
      lvl_opts[[i]] <- list(
        level = i,
        dataLabels = list(enabled = FALSE))
    }
  }
  
  # Create treemap object, to save as png and html later (outside this func)
  tm <- summ %>%
    data_to_hierarchical(c(...), n, brewer.pal(n = 8, name = "Dark2")) %>%
    hchart(
      type = "treemap",
      allowTraversingTree = TRUE,
      levelIsConstant = FALSE,
      levels = lvl_opts,
      drillUpButton = list(text = "←")) %>%
    hc_chart(
      style = list(fontFamily = "Lexend")) %>%
    hc_tooltip(
      pointFormat = "<b>{point.name}</b>: {point.value} samples<br/>",
      useHTML = TRUE) %>%
    hc_title(
      text = tm_title,
      align = "left") %>%
    hc_subtitle(
      text = tm_subtitle,
      align = "justify") %>% 
    hc_caption(
      text = tm_caption,
      align = "justify") %>%
    hc_credits(
      enabled = TRUE, text = "Data Source: GISAID (2020-2023)",
      href = "https://gisaid.org/")
  
  return(tm)
}

# Workaround for broken saveWidget
saveWidget2 <- function(widget, file) {
  htmlwidgets::saveWidget(widget, file = "tmp.html", selfcontained = TRUE)
  file.rename("tmp.html", file)
}

# Confirm duplicates in genome samples
confirm_fasta_dupes <- function(fasta_all) {
  for (i in which(duplicated(fasta_all))) {
    test <- fasta_all[i]
    for (j in which(duplicated(fasta_all))) {
      if (i == j) next
      else if (isTRUE(ape::all.equal.DNAbin(fasta_all[[i]], fasta_all[[j]])))
        print(sprintf("FOUND! fasta_all[%d] (%s) matches fasta_all[%d] (%s)", i, names(fasta_all[i]), j, names(fasta_all[j])))
    }
  }
}

# Get the lineages in each variant and write to a csv file and txt files
get_var_lin <- function(metadata_all, write_path = "data/overview") {
  df <- metadata_all %>%
    dplyr::select(pangolin_lineage, variant) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(pangolin_lineage)
  write_path <- paste0(write_path, "/vars")
  if (!dir.exists(write_path))
    dir.create(write_path)
  write_csv(df, paste0(write_path, "/var_lin.csv"))
  for (var in unique(metadata_all$variant)) {
    output <- df %>%
      dplyr::filter(variant == var) %>%
      dplyr::select(pangolin_lineage) %>%
      dplyr::arrange()
    write_lines(unlist(output),
                paste0(write_path, sprintf("/%s_lineages.txt",
                        stringr::str_to_lower(
                          stringr::str_replace_all(var, " ", "")))))
  }
}

# experimental plotly treemap
treemap2 <- function(metadata_all) {
  df <- metadata_all %>%
    dplyr::mutate(dplyr::across(c(variant, division_code, pangolin_lineage, strain), as.character), .keep="used")
  
  data <- df %>%
    dplyr::mutate(root = "Region") %>%
    dplyr::select(labels = division_code, parents = root) %>%
    dplyr::add_row(labels = df$variant, parents = df$division_code) %>%
    dplyr::add_row(labels = df$pangolin_lineage, parents = df$variant) %>%
    dplyr::add_row(labels = df$strain, parents = df$pangolin_lineage) %>%
    dplyr::distinct()
  
  fig <- plot_ly(
    type="treemap",
    labels=data$labels,
    parents=data$parents,
    values=rep(1,times=nrow(df))
  )
  fig
  
}
# File: helper.R
# To contain auxiliary functions for all code/R/ files.

# Appends system information and parameters used to a text file.
# specified by output_path and filename.
# This function only supports Windows 10/11 and Linux.
# In Windows, assumes that there's only one processor installed.
write_to_log <- function(output_dir, filename, log_string) {
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
    specs <- c(systeminfo, cpuinfo)
    write_lines(c(as.character(Sys.time()), specs, log_string,
                  "------\n"), output_path, append = TRUE)
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
    write_lines(c(as.character(Sys.time()), specs, log_string,
                  "------\n"), output_path, append = TRUE)
  } else {
    write_lines(c(as.character(Sys.time()),
                  "OS not supported by logger!", log_string,
                  "------\n"), output_path, append=TRUE)
  }
  close(fileConn)
}

# Returned time is of the format YYYYMMDDHHMMSSNNN
get_time <- function() {
  ret <- gsub("[.:-]|\\s", "", as.character(Sys.time()))
  ret <- substring(ret, 1, nchar(ret)-3)
}

# MAKE TREEMAPS
# Plot treemaps with appropriate drilldowns using highcharter.
make_treemaps <- function(metadata_all, write_path, stamp) {
  # Wrapper function for making treemaps with ...length() levels
  treemap <- function(df, ..., write_path, stamp) {
    # Create df containing the columns to summarize.
    summ <- df %>% select(...) %>% tibble(n = rep(1, nrow(df)))
    # Generate level JSONs
    lvl_json <- list()
    for (i in 1:...length()) {
      if (i == 1) {
        lvl_json[[i]] <- list(
          level = 1,
          dataLabels = list(enabled = TRUE,
                            format = "{point.name}<br>{point.value}"),
          borderColor = "white",
          borderWidth = 1)
      } else if (i == 2) {
        lvl_json[[i]] <- list(
          level = 2,
          dataLabels = list(enabled = TRUE,
                            style = list(fontSize = "0.8em"))
        )
      } else {
        lvl_json[[i]] <- list(
          level = i,
          dataLabels = list(enabled = FALSE)
        )
      }
    }
    
    # Create treemap object, to save as png and html later (outside this func)
    tm <- summ %>%
      data_to_hierarchical(c(...), n) %>%
      hchart(type = "treemap",
             allowTraversingTree = TRUE,
             levelIsConstant = FALSE,
             levels = lvl_json)
    tm
  }
  
  metadata_all %>% treemap(variant, ph_region, pangolin_lineage)
}


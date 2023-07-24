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

# Plot treemap with appropriate drilldowns using highcharter
make_treemap <- function(metadata) {
  set.seed(1234)
  
  # Summary table
  summary.table <- metadata_all %>% 
    group_by(variant) %>% 
    summarise(
      nb_variant = n(), 
      nb_division_exposure = length(unique(division_exposure))
    ) %>% 
    arrange(-nb_variant, -nb_division_exposure)
  summary.table
  
  hc <- summary.table %>%
    hchart(
      "treemap", 
      hcaes(x = variant, value = nb_variant)
    ) %>%
    hc_add_theme(hc_theme_ggplot2()) %>%
    hc_colorAxis(stops = color_stops(colors = viridis::inferno(n = 6, direction = -1)),
                 max =  50)
  hc
  
  hc <- data_to_hierarchical(metadata_all, c(division_exposure, variant), division_exposure)
  hchart(hc, type = "treemap")
  
  set.seed(110)
  
  ex <- data.frame(
    l1 = metadata_all$division_exposure,
    l2 = metadata_all$variant,
    l3 = metadata_all$pangolin_lineage,
    count = rep(1,nrow(metadata_all))
  )
  
  ex %>% 
    data_to_hierarchical(c(l1, l2, l3), count) %>%
    hchart(type = "treemap",
           allowTraversingTree = TRUE,
           levelIsConstant = FALSE,
           levels = list(
             list(level = 1, dataLabels = list(enabled = TRUE,
                                               format = "{point.name}<br>
                                               {point.value}"),
                  borderColor = "white", borderWidth = 1),
             list(level = 2, dataLabels = list(enabled = TRUE,
                                               style = list(fontSize = "0.8em"))),
             list(level = 3, dataLabels = list(enabled = FALSE))
           )
    )
}


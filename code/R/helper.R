# File: helper.R
# To contain auxiliary functions for all code/R/ files.

# Appends system information and parameters used to a text file.
# specified by output_path and filename.
# This function only supports Windows 10/11 and Linux.
# Assumes that there's only one processor installed.
write_to_log <- function(output_dir, filename, log_string) {
  if (!dir.exists(output_dir))
    dir.create(output_dir)
  output_path <- paste(output_dir, filename, sep = "/")
  
  fileConn<-file(output_path)
  
  if (pacman::p_detectOS() == "Windows") {
    systeminfo <- system("systeminfo", intern=TRUE)
    systeminfo <- systeminfo[c(3,4,13,14,15,16,17,25)]
    cpuinfo <- system("WMIC CPU Get DeviceID,NumberOfCores,NumberOfLogicalProcessors",
                      intern = TRUE)
    cpuinfo <- cpuinfo[-length(cpuinfo)]
    specs <- c(systeminfo, cpuinfo)
    write_lines(c(as.character(Sys.time()), specs, log_string,
                  "------\n"), output_path, append=TRUE)
  } else if (pacman::p_detectOS() == "Linux" && 
             (Sys.getenv("DOCKER_RUNNING") == "" ||
              Sys.getenv("DOCKER_RUNNING") == FALSE)) {
    device <- paste(system(paste("cat /sys/devices/virtual/dmi/id/sys_vendor",
                                 "/sys/devices/virtual/dmi/id/product_name",
                                 "/sys/devices/virtual/dmi/id/product_version"),
                           intern=T), collapse=' ')
    lsb_release <- system("lsb_release -a", intern=T)[c(3,5)]
    lscpu <- system("lscpu", intern=T)[c(1,5,8,11,16,17,30,33,34,35)]
    mem <- system("grep MemTotal /proc/meminfo", intern=T)
    specs <- c(device, lsb_release, lscpu, mem)
    write_lines(c(as.character(Sys.time()), specs, log_string,
                  "------\n"), output_path, append=TRUE)
    # If DOCKER_RUNNING is True, we're automatically on Linux Container
  } else if (pacman::p_detectOS() == "Linux") {
  } else {
    write_lines(c(as.character(Sys.time()),
                  "OS not supported by logger!", log_string,
                  "------\n"), output_path, append=TRUE)
  }
  close(fileConn)
}

# Returned time is a file function
get_time <- function() {
  ret <- gsub("[.:-]|\\s", "", as.character(Sys.time()))
  ret <- substring(ret, 1, nchar(ret)-3)
}
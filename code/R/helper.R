# File: helper.R
# Supports: kmer-analysis.R, preprocess.R
# To contain auxiliary functions

if (!require("pacman"))
  install.packages("pacman", repos = "https://cran.case.edu")
library(pacman)

# Write system information and parameters to a text file specified in
# output_path. Creates text file in output_path if not exists.
# Overwrites if exists. Only works for Windows and Linux.
# Assumes that there's only one processor installed.
paramsLog <- function(outputDir, filename, paramString) {
  if (!dir.exists(outputDir)) dir.create(outputDir)
  output_path <- paste(outputDir, filename, sep = '/')
  fileConn<-file(output_path)
  if (pacman::p_detectOS() == 'Windows') {
    specs <- system('systeminfo', intern=TRUE)[
      c(3,4,8,13,14,15,16,17,length(specs))]
    # WMIC CPU Get DeviceID,NumberOfCores,NumberOfLogicalProcessors
    # Get-WmiObject win32_VideoController
    write_lines(c(as.character(Sys.time()), specs, paramString,
                  '------\n'), output_path, append=TRUE)
  } else {
    hostnamectl <- trimws(system('hostnamectl', intern=T)[c(6,7,9,10)])
    lsb_release <- system('lsb_release -a', intern=T)
    lsb_release <- gsub('\\t', ' ', lsb_release[length(lsb_release)])
    lscpu <- system('lscpu', intern=T)[c(1,5,8,11,16,17,30,33,34,35)]
    mem <- system('grep MemTotal /proc/meminfo', intern=T)
    specs <- c(hostnamectl, lsb_release, lscpu, mem)
    write_lines(c(as.character(Sys.time()), specs, paramString,
                  '------\n'), output_path, append=TRUE)
  } 
  close(fileConn)
}

timeString <- function() {
  ret <- gsub('[.:-]|\\s', '', as.character(Sys.time()))
  ret <- substring(ret, 1, nchar(ret)-3)
  return(ret)
}
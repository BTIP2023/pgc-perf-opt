# File: helper.R
# Supports: kmer-analysis.R, preprocess.R
# To contain auxiliary functions

# Write system information and parameters to a text file specified in
# output_path. Creates text file in output_path if not exists.
# Overwrites if exists. Only works in Windows.
paramsLog <- function(output_path, paramString) {
  fileConn<-file(output_path)
  specs <- system('systeminfo', intern=TRUE)
  write_lines(c(as.character(Sys.time()),
               specs[3], specs[4], specs[8],
               specs[13], specs[14], specs[15],
               specs[16], specs[17], specs[25],
               paramString,
               '------\n'), output_path, append=TRUE)
  close(fileConn)
}

timeString <- function() {
  ret <- gsub('[.:-]|\\s', '', as.character(Sys.time()))
  ret <- substring(ret, 1, nchar(ret)-3)
  return(ret)
}
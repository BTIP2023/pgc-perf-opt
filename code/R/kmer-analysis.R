# File: kmer-analysis.R
# Main function: get_kmers
# Generates kmer_k_stamp.csv files from fasta and meta using kcount.
# Generated files will have the metadata cols appended, starting from strain.

get_kmers <- function(fasta, metaData, k, stamp) {
  message(sprintf("\nPerforming %d-mer analysis... ", k), appendLF = FALSE)
  kmers <- kcount(fasta, k = k)
  kmer_df <- data.frame(kmers)
  message("DONE.")
  
  # Append metaData
  kmer_df <- cbind(kmer_df, metaData)
  
  # Write to a csv file in data/kmers/
  if (!dir.exists("data/kmers"))
    dir.create("data/kmers")
  output_dir <- paste("data/kmers", sprintf("kmer_%d_%s.csv", k, stamp),
                      sep = "/")
  
  message(paste0("Writing kmer data to ", output_dir, "... "), appendLF = FALSE)
  write_csv(kmer_df, output_dir)
  message(paste0("Done writing to", output_dir))
  
  message("Writing logs... ", appendLF = FALSE)
  write_to_log(output_dir = "data/kmers", filename = "log.txt",
               log_string = sprintf("timestamp = %s\nseed = %d, strat_size = %d, k-value = %d",
                                    stamp, seed, strat_size, k))
  message("Done writing logs.")
}
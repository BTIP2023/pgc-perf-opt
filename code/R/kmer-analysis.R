# File: kmer-analysis.R
# Main function: get_kmers
# Generates kmer_k_stamp.csv files from fasta and meta using kcount.
# Generated files will have the metadata cols appended, starting from strain.

get_kmers <- function(fasta, metaData, k, stamp) {
  print(sprintf("Performing %d-mer analysis...", k))
  
  kmers <- kcount(fasta, k = k)
  kmer_df <- data.frame(kmers)
  
  # Append metaData
  kmer_df <- cbind(kmer_df, metaData)
  
  # Write to a csv file in data/kmers/
  if (!dir.exists("data/kmers"))
    dir.create("data/kmers")
  output_dir <- paste("data/kmers", sprintf("kmer_%d_%s.csv", k, stamp),
                      sep = "/")
  
  print(paste("Writing kmer data to", output_dir))
  write.csv(kmer_df, output_dir)
  
  write_to_log(output_dir = "data/kmers", filename = "log.txt",
               log_string = sprintf("timestamp = %s\nseed = %d, strat_size = %d, k-value = %d",
                                    stamp, seed, strat_size, k))
}
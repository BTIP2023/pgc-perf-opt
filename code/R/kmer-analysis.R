# File: kmer-analysis.R
# Main function: get_kmers
# Generates kmer_k_stamp.csv files from fasta and meta using kcount.
# Generated files will have the metadata cols appended, starting from strain.

get_kmers <- function(fasta, metaData, k, stamp) {
  message(sprintf("Performing %d-mer analysis... ", k), appendLF = FALSE)
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
  message(paste0("Writing kmer data to ", output_dir, "... DONE."))
}

generate_kmer_heatmap <- function(kmers, results_path, k) {
  # Process kmers dataframe
  df <- kmers
  
  # Preprocessing of Data for heatmap generation
  print("Preprocessing data for heatmap generation...")
  slice_col <- which(colnames(df) == "strain")
  df <- df[, 2:(slice_col)]
  
  df_long <- pivot_longer(data = df, 
                          cols = -c(strain), 
                          names_to = "kmer",
                          values_to = "frequency") %>%
    mutate(freq_norm = (frequency-min(frequency))/(max(frequency)-min(frequency)))
  
  # Heatmap generation
  print("Generating kmer heatmap...")
  p <- ggplot (df_long, aes(x = kmer, y= strain, fill = freq_norm)) + 
    geom_tile() +
    xlab(label = "kmer") +
    ylab(label = "Sample") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    scale_fill_gradient (low = "#00AFBB", high = "#FC4E07")
  
  # Saving heatmap as PNG, HTML, and RData
  results_path <- paste0(results_path, "/heatmaps")
  save_plot("heatmap", results_path, k, p)
}

# k-mer descriptors -> n-gram descriptors, use k = 7 for final output
# For each sample, then for all samples (mean)
# I'm making this in shiny because my gosh there are 23000 plus samples
# It would be horrifying to check out all of those manually
generate_kmer_wordcloud <- function(kmer_df, write_path = "data/kmers", k) {
  write_path <- paste0(write_path, "/wordclouds")
  if(!dir.exists(write_path))
    dir.create(write_path)
  
}
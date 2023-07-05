
#Import all the required libraries
library(ape)
library(kmer)
library("dplyr")

#set the working directory 
setwd("")

#Function computes the kmer of given length 
kmer_df = function(filePath,variant,k){
  kmer = read.FASTA(filePath)
  kmer3 = kcount(kmer, k = k)
  target = rep((variant),length(kmer))
  kmer_df1 = data.frame(kmer3, target)
  return(kmer_df1)
}

#Different lengths of kmers to be used 
kmer_list = list(3,5,7)

for(k in kmer_list)
{
  alpha = kmer_df('alpha.fasta','alpha',k)
  beta = kmer_df('beta.fasta','beta',k)
  gamma = kmer_df('gamma.fasta','gamma',k)
  omicron = kmer_df('omicron.fasta','omicron',k)
  delta = kmer_df('delta.fasta','delta',k)
  
  data = bind_rows(alpha, beta, delta, gamma, omicron)
  outputName = sprintf("covid_kmer_%d.csv",k)
  #store the combined data in a csv file 
  write.csv(data, outputName)
}
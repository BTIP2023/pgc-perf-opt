# Code Guide

## For kmer-analysis.R
* `seed` variable is there for reproducibility of random values. You may set this.
* You may also set `sampleSize` to desired number of random samples.
* This code outputs the generated kmers with metadata to `data/kmers`. 
* Note: warnings that appear (parsing errors due to conflicting column types and values, and NAs due to type coercion) can be ignored as these are handled in the code.
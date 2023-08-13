# Data
`data/GISAID` contains the raw GISAID data to be fed to the pipeline's entry point, `code/R/preprocess.R`.
If `write_fastacsv <- TRUE` in `code/pipeline.R` files, `data/interm` will additionally contain intermediate but thoroughly sanitized FASTA and metadata files of the entire[^1] and the subsampled dataset.

As such, `data/interm` can be used to fix value issues in the raw GISAID datasets.
See `get_sample` and `sanitize_sample` in code/R/preprocess.R for all the relevant issues and fixes.
Moreover, `data/overview` contains overviews (in .txt, .csv, and plots) of the entire[^1] and the subsampled dataset.

## On Timestamps
Generated files in `data/interm`, `data/overview`, and `data/kmers` are timestamped according to
when they were generated. The timestamps have the format `YYYYMMDDHHMMSSNNN`. Note that
a timestamp is generated at the top of the pipeline and only passed to succeeding functions.
Some pipeline functions can also work without getting a timestamp by ordering the files
and selecting the most recent one.

Untimestamped files[^1] can also be used (i.e. files suffixed by stratum sizes) by setting `stamp <- str(strat_size)` where `str(strat_size)` is the suffix of the desired untimestamped file to use (which may or may not be pre-generated). Just make sure to also set `strat_size` to whatever you're setting `stamp` to!

[^1]: Untimestamped files in `data/interm`, `data/overview`, and `data/kmers` are
files that contain generated data pertaining to the ENTIRE dataset (all samples in GISAID included).

## On Sharing and Using GISAID Data
For data obtained from GISAID, only **accession numbers** will be included in this and the root README files' references. GISAID (2012) gives the following reason:

> GISAID does not promote the release of data to databases where access to data is anonymous and the rights of the submitter are relinquished.  GISAID already provides the public with open access to data in a transparent way.

All GISAID data are stored in `data/GISAID`, but the contents won't be committed (via `.gitignore`).

### GISAID Directory Usage
To main devs: To use the GISAID directory, open the Google Drive folder shared by the mentor, and then navigate to `data/GISAID`. Copy the contents of that (`.tar` files) into your local repo's `data/GISAID`. The aforementioned directory has to be updated every time there's a change in the data.

To future researchers: With the complete list of Accession Numbers at hand (see [`data/overview/accession.txt`](`data/overview/accession.txt`)), and with your own GISAID Access Credentials, you may be able to re-download the GISAID data used for our analysis from [gisaid.org](https://gisaid.org/).

---
## References
GISAID. (2012). FAQ. https://gisaid.org/help/faq/

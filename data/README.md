# Data
This directory will contain the raw `.fasta` files to be processed using k-mer analysis, generating `.csv` files.

## On GISAID Data
For data obtained from GISAID, only **accession numbers** will be included in this and the root README files' references. GISAID (2012) gives the following reason:

> GISAID does not promote the release of data to databases where access to data is anonymous and the rights of the submitter are relinquished.  GISAID already provides the public with open access to data in a transparent way.

All GISAID data are stored in `~/data/GISAID`, but the contents won't be committed (via `.gitignore`).

### GISAID Directory Usage
To devs: To use the GISAID directory, open the Google Drive folder shared by the mentor, and then navigate to data/GISAID/datasets. Copy the contents of that (`.fasta` and `.tsv`) into your local repo's `~/data/GISAID/sequences`. The aforementioned directory has to be updated every time there's a change in the data.

---
## References
GISAID. (2012). FAQ. https://gisaid.org/help/faq/

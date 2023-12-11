# meth

This project transfers raw methylation illumina / 450K data taken from here
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134089

Into beta values using minfi
https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html

And then calculating DNAmAge using multi clocks from here
https://www.bioconductor.org/packages/devel/bioc/vignettes/methylclock/inst/doc/methylclock.html

------------------- 12.11.23 -

https://github.com/mkholod/meth/blob/main/DNAmAge_calc_full.R
* Reads the meth. arrays and converts them to beta.csv file
* Calculates myDNAmAge_with_acceleration_age

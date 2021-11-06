###########################################################################################
# ChAmp https://bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
# alternatively - but I didn't run the one above was enough
# BiocManager::install(c("minfi","ChAMPdata","Illumina450ProbeVariants.db","sva","IlluminaHumanMethylation450kmanifest","limma","RPMM","DNAcopy","preprocessCore","impute","marray","wateRmelon","goseq","plyr","GenomicRanges","RefFreeEWAS","qvalue","isva","doParallel","bumphunter","quadprog","shiny","shinythemes","plotly","RColorBrewer","DMRcate","dendextend","IlluminaHumanMethylationEPICmanifest","FEM","matrixStats","missMethyl","combinat"))

library("ChAMP")

# this is a test dir, once done testDir == "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/ChAMPdata/extdata"
testDir=system.file("extdata",package="ChAMPdata") 
myLoad <- champ.load(testDir, arraytype="450K") # myLoad$beta has the beta values

##############
# now our example - clear previous values
rm(testDir) 
rm(myLoad)
rm(probe.features)
rm(multi.hit)
rm(hm450.manifest.hg19)
rm(Anno)

testDir <- "./samples"
#myLoad <- champ.load(testDir, arraytype = "450K") # ok that failed on filtering so replacing with champ.import(testDir) + champ.filter()
myImport <- champ.import(testDir)

# failing
# Error in `[<-`(`*tmp*`, type_II, i, value = numeric(0)) : 
#  no 'dimnames' attribute for array

tirosh_idat <- illuminaio::readIDAT("./samples/GSM3936886_202163530086_R08C01_Red.idat")
champ_idat_sample <- illuminaio::readIDAT("/Library/Frameworks/R.framework/Versions/4.1/Resources/library/ChAMPdata/extdata/7990895118_R01C01_Red.idat")


source("./temp/champ_import.r")
misha_champ_import(testDir)

# line 125, M is false for all
# M <- IDAT[Anno$Annotation[, "M.index"][M.check], ]
# message("    Extracting UnMeth Matrix...")

################# didn't go beyond this line #############
myLoad <- champ.filter() # that worked
# CpG.GUI() # opens a window - doesn't show anything
champ.QC() # didn't work
champ.process(directory = testDir)
# myLoad <- champ.load(testDir)



###########
g_new_rownames = c("111","222")
misha <- sapply(g_new_rownames,function(x) gsub('^([0-9]{1})([0-9]+)$', '\\10\\2', as.character(x)))

############
## check values in Anno - when sorted the start with [1] "G-10600322" and our samples start with 1600303 we need to add 9,000,000
m_index <- Anno$Annotation[, "M.index"]
sorted_m_index <- sort(m_index)


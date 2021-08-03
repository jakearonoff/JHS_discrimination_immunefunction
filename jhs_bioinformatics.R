library(meffil)
library(tidyverse)
library(FlowSorted.Blood.EPIC) # uses minfi 
library(IlluminaHumanMethylationEPICmanifest)
library(ewastools)
library(stringi)
library(data.table)
library(svd)

### NOTE: this analysis required a very large amount of computer memory, and as a result the code for estimating 
#         lymphocyte %'s and the betas for later estimation of aging indicators was run through 11 batches. 
#         Normally analyzing methylation in batches is problematic due to potential systematic batch differences. 
#         However, the methods used here for estimating cell type %'s and aging indicators are robust to potential 
#         batch effects. 

###################################
# Contamination Check using ewastools package (this was run on the whole sample at once, not in batches)
###################################

targets_qc = read.metharray.sheet("path_to_directory_with_sample_sheets_and_IDATS")
pheno = targets_qc # storing sample sheet info 
meth = read_idats(pheno$Basename,quiet = FALSE) # reading in methylation data from IDATS
beta = meth %>% detectionP %>% mask(0.01) %>% correct_dye_bias %>% dont_normalize
snps = meth$manifest[probe_type=="rs",index]
snps = beta[snps,]
genotypes = call_genotypes(snps,learn=FALSE)
pheno$outlier = snp_outliers(genotypes)
contam <- pheno[pheno$outlier > -4,c("Sample_Name")]
contam # potentially contaminated samples that were removed prior to quality control procedures 


###################################
# QC using meffil package (this was run on the whole sample at once, not in batches)
###################################
targets_qc = read.metharray.sheet("path_to_directory_with_sample_sheets_and_IDATS")
targets_qc <- filter(targets_qc, !(Sample_Name %in% contam)) # removing potentially contaminated samples
targets_qc$Sex <- ifelse(targets_qc$Gender == "Male","M","F") # recoding sex variable for QC pipeline 
qc.objects <- meffil.qc(targets_qc)
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1, 
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5
) # defining parameters for QC
qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters
) # performing QC
meffil.qc.report(qc.summary, output.file="filepath.html")# output QC results into html
badsamples <- as.data.frame(qc.summary$bad.samples) 
badsamplesnames <- badsamples$sample.name # store names of samples that didn't pass QC


###################################
# Estimating % of different lymphocyte types present in PBMC's using FlowSorted.Blood.EPIC package 
###################################

data(IDOLOptimizedCpGs) # get IDOL CpG's (probes used for deconvolution)
targets = read.metharray.sheet("path_to_directory_with_sample_sheets_and_IDATS")
targets <- filter(targets, !(Sample_Name %in% contam)) # removing potentially contaminated samples
targets <- filter(targets, !(Sample_Name %in% badsamplesnames)) # remove samples that failed QC 
rgSet <- read.metharray.exp(targets = targets, force = TRUE) # reading in methylation data from IDAT files 
wbc <- estimateCellCounts2(rgSet, compositeCellType = "Blood",
                           processMethod = "preprocessNoob",
                           probeSelect = "IDOL",
                           cellTypes = c("CD8T", "CD4T", "Bcell"),
                           referencePlatform =
                             "IlluminaHumanMethylationEPIC",
                           referenceset = NULL,
                           IDOLOptimizedCpGs =IDOLOptimizedCpGs,
                           returnAll = TRUE)    #  estimates cell type %'s
write.csv(wbc$counts, file = "filepath.csv") # write results to a .csv


###################################
# Outputting beta values to a .csv file to input later into Horvath and colleagues website for aging estimates 
#  website url: http://dnamage.genetics.ucla.edu   (used new clock)
###################################

rcobj <- ratioConvert(wbc$normalizedData, what = "beta", keepCN = F)
noob.betas <- getBeta(rcobj)
class(noob.betas)
noob.betasd <- as.data.frame(noob.betas)
str(noob.betasd)
noob.betasd$Name <- rownames(noob.betasd) # Make rownames (which are probe names) a column in the matrix
probes <- read.csv("path/datMiniAnnotation3.csv") # read in Horvath's probe list, downloaded from website 
fil.noob.betas <- merge(probes, noob.betasd, by = "Name", all.x = TRUE, all.y = FALSE) # get betas for Horvath's probe list
nrow(fil.noob.betas) # should be 30084 
fil.noob.betas <- fil.noob.betas[-c(2:7)] # removing extra columns
colnames(fil.noob.betas)[1] <- "ProbeID"
write.csv(fil.noob.betas,"filepath.csv") # save sample probes into .csv for later upload to Horvath website 












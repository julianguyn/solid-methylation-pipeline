# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(limma)
    library(ggplot2)
    library(ggpubr)
    library(RColorBrewer)
    library(viridis)
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(rGREAT)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

mat <- readRDS("data/rawdata/cholangio/MeDIP_Solid_cholangiocarcinoma_log2CPM.rds")

pheno <- fread("data/rawdata/cholangio/Pheno_cholangio_all_IDHwt_Vs_IDHmut.txt") #64
#pheno2 <- fread("data/rawdata/cholangio/Pheno_cholangio_baseline_IDHwt_Vs_IDHmut.txt") #16

# keep only cholangio samples
mat <- mat[,colnames(mat) %in% pheno$MeDIP_ID]
pheno <- pheno[match(colnames(mat), pheno$MeDIP_ID),] # 38

ids <- unique(pheno$solid_patient)
for (id in ids) {
    subset <- pheno$response[pheno$solid_patient == id]
    cat(id, ":", unique(subset), "\n")
}

###########################################################
# Quick PCA
###########################################################

# rank date for each sample
pheno$date <- as.Date(pheno$blood_collection_date, format = "%Y %b %d")
pheno <- pheno[order(pheno$solid_patient, pheno$date), ]
pheno$timepoint <- unlist(lapply(as.numeric(table(pheno$solid_patient)), seq_len))
pheno$timepoint <- factor(paste0("T", pheno$timepoint), levels = paste0("T", 1:12))

# add final response
baseline <- pheno[pheno$response == "baseline",]
baseline$final_response <- ifelse(
    baseline$solid_patient %in% pheno$solid_patient[pheno$response == "PD"],
    "PD", "SD"
)
# pheno[,c(2,6, 10, 17)]
# baseline[,c(2, 6, 10, 17, 18)]
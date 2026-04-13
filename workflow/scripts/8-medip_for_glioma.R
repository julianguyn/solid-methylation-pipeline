# quickl subset medip for glioma
# messy code!

library(data.table)


mat <- readRDS("data/rawdata/cholangio/MeDIP_Solid_cholangiocarcinoma_log2CPM.rds")
pheno <- fread("data/rawdata/cholangio/Pheno_cholangio_all_IDHwt_Vs_IDHmut.txt") #64

# load in all medip
medip <- readRDS("data/rawdata/MeDIP_Solid_all_log2CPM.rds")

# keep only cholangio samples
cholangio <- intersect(pheno$MeDIP_ID, colnames(mat))

# remove cholangio
pheno <- pheno[-which(pheno$MeDIP_ID %in% common),]

# get glioma samples
load("data/procdata/RGSets_targets.RData")

need <- targets$Subject
need <- paste0("Solid_", sub(".*-", "", need))
need <- gsub("_0", "_", gsub("_0", "_", need))

have <- paste0("Solid_", sub("batch.*_", "", colnames(medip)))
have <- sub("-.*", "", have)
have <- have[order(have)]

# gliomas with medip
glioma <- need[need %in% have] #20
missing <- need[-which(need %in% have)]

## TODO: subset matrix
medip <- medip[,colnames(medip) %in% pheno$MeDIP_ID]

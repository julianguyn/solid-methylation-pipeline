# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(limma)
    library(ggplot2)
    library(ggpubr)
    library(RColorBrewer)
    library(viridis)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

mat <- readRDS("data/rawdata/cholangio/MeDIP_Solid_cholangiocarcinoma_log2CPM.rds")

pheno1 <- fread("data/rawdata/cholangio/Pheno_cholangio_all_IDHwt_Vs_IDHmut.txt") #64
pheno2 <- fread("data/rawdata/cholangio/Pheno_cholangio_baseline_IDHwt_Vs_IDHmut.txt") #16

# keep only PD & SD
pheno <- pheno1[pheno1$response != "baseline",] # 38

###########################################################
# Quick PCA
###########################################################

# rank date for each sample
pheno$date <- as.Date(pheno$blood_collection_date, format = "%Y %b %d")
pheno <- pheno[order(pheno$solid_patient, pheno$date), ]
pheno$timepoint <- unlist(lapply(as.numeric(table(pheno$solid_patient)), seq_len))
pheno$timepoint <- factor(paste0("T", pheno$timepoint), levels = paste0("T", 1:12))

# perform PCA
mat <- mat[,match(pheno$MeDIP_ID, colnames(mat))] |> as.matrix()
pca_res <- prcomp(t(mat))$x

p1 <- ggplot(pca_res, aes(x = PC1, y = PC2, color = pheno$solid_patient, shape = pheno$response)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "Set3") +
    theme_bw()

p2 <- ggplot(pca_res, aes(x = PC1, y = PC2, color = pheno$timepoint, shape = pheno$response)) +
    geom_point(size = 3) +
    scale_color_viridis_d(option = "rocket", begin = 0.2) +
    theme_bw()

p <- ggarrange(p1, p2, ncol=2)
ggsave("data/results/figures/cholangio/pca_dmrs_medip.png", p, width = 11.5, height = 4.5)

###########################################################
# Run linear model (SD-PD)
###########################################################

# set variables
labels <- factor(pheno$response)

# model
design <- model.matrix(~ labels) # 0 DMRs with FDR<0.05
fit <- lmFit(mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1)

###########################################################
# Run linear model (patient + SD-PD)
###########################################################

# set variables
labels <- factor(pheno$response)
patientid <- factor(pheno$solid_patient)

# model
design <- model.matrix(~ patientid + labels) # 1 DMR with FDR<0.05
fit <- lmFit(mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1)

###########################################################
# Run linear model (corfit)
###########################################################

# design matrix
design <- model.matrix(~ labels)

# estimate correlation within patients
corfit <- duplicateCorrelation(mat, design, block = patientid)

# model (0 DMRs)
fit <- lmFit(mat, design, block = patientid, correlation = corfit$consensus)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1)

###########################################################
# Keep only first SD and first PD
###########################################################

table(pheno$response) # PD: 21, SD: 17

# keep first SDs
sd_pheno <- pheno[pheno$response == "SD",]
idx <- cumsum(as.numeric(table(sd_pheno$solid_patient)))
idx <- c(1, idx+1)[-c(length(idx)+1)]
sd_pheno <- sd_pheno[idx,]

# keep first PDs
pd_pheno <- pheno[pheno$response == "PD",]
idx <- cumsum(as.numeric(table(pd_pheno$solid_patient)))
idx <- c(1, idx+1)[-c(length(idx)+1)]
pd_pheno <- pd_pheno[idx,]

###########################################################
# Run linear model (SD-PD)
###########################################################

# subset matrix
subset_pheno <- rbind(sd_pheno, pd_pheno)
subset_mat <- mat[,match(subset_pheno$MeDIP_ID, colnames(mat))]

# set variables
labels <- factor(subset_pheno$response)

# model
design <- model.matrix(~ labels) # 0 DMRs with FDR<0.05
fit <- lmFit(subset_mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1) # 17 DMRs

###########################################################
# Run linear model (patient + SD-PD)
###########################################################

# set variables
labels <- factor(subset_pheno$response)
patientid <- factor(subset_pheno$solid_patient)

# model
design <- model.matrix(~ patientid + labels) # 1 DMR with FDR<0.05
fit <- lmFit(subset_mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1) #479 DMRs

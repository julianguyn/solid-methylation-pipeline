# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(minfi)
    library(limma)
    library(DMRcate)
    library(ggplot2)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

message("Loading in data")

# load in RGSets and targets
load("data/procdata/RGSets_targets.RData")
targets$match <- gsub(".*/", "", targets$Basename)

# get beta values
bVals <- fread("data/results/data/bVals.tsv", data.table = FALSE)
rownames(bVals) <- bVals$V1
bVals$V1 <- NULL

# get m values
mVals <- fread("data/results/data/mVals.tsv", data.table = FALSE)
rownames(mVals) <- mVals$V1
mVals$V1 <- NULL

# get annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- anno[rownames(mVals), ]

###########################################################
# Format targets
###########################################################

# quick sanity check
table(colnames(mVals) == colnames(bVals))
targets <- targets[match(colnames(mVals), targets$Subject),]
table(colnames(mVals) == targets$Subject)

targets$responder <- ifelse(
    targets$Response_RANO %in% c("CR", "PR"),
    "Responder",
    "NonResponder"
)
targets$responder <- factor(targets$responder)
targets$Sex <- factor(targets$Sex)


###########################################################
# DMP analysis
###########################################################

# no DMPs at adj.P.Val < 0.05 for 4 responders vs non-responders

# design matrix and model
design <- model.matrix(~ responder + Age + Sex, data = targets)
fit <- lmFit(mVals, design)
fit2 <- eBayes(fit)

# get dmps
dmps <- topTable(fit2, coef = 2, number = Inf, adjust.method = "BH", sort.by = "p")

###########################################################
# DMR analysis
###########################################################

# get cpg annotation
myannotation <- cpg.annotate(
    datatype = "array",
    object = as.matrix(mVals),
    what = "M",
    analysis.type = "differential",
    design = design,
    contrasts = FALSE,
    coef = "responderResponder",
    arraytype = "EPICv1",
    fdr = 0.1
)

# get dmrs
dmrs <- extractRanges(
    dmrcate(myannotation, lambda = 1000, C = 2), genome = "hg19"
)
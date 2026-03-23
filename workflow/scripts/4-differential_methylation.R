# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(minfi)
    library(limma)
    library(DMRcate)
    library(ggplot2)
    library(patchwork)
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

# load in updated metadata
meta <- read.csv("metadata/clin_updated.csv")

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

targets <- targets[match(colnames(mVals), targets$Subject),]
meta <- meta[match(colnames(mVals), meta$Subject),]

targets$responder <- ifelse(
    targets$Response_RANO %in% c("CR", "PR"),
    "Responder",
    "NonResponder"
)
targets$responder <- factor(targets$responder)
targets$Sex <- factor(targets$Sex)
targets$CDKN2 <- factor(meta$CDKN2A_B_mut)

###########################################################
# Target evaluation
###########################################################

toPlot <- targets[order(targets$PFS_Days, decreasing = TRUE),]
toPlot$Subject <- factor(toPlot$Subject, levels = toPlot$Subject)
toPlot$Response_RANO <- factor(toPlot$Response_RANO, levels = names(mrecist_pal))

p1 <- ggplot(toPlot, aes(x = Subject, y = 1, fill = CDKN2)) +
    geom_tile(color = "black") +
    scale_fill_manual(values = c("#2F4858", "#86BBD8"), na.value = "white") +
    theme_void()

p2 <- ggplot(toPlot, aes(x = Subject, y = PFS_Days, fill = Response_RANO)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("RECIST", values = mrecist_pal) +
    theme_minimal() +
    theme(
        panel.border = element_rect(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank()
    ) +
    labs(y = "PFS (Days)")

filename <- "data/results/figures/targets/targets.png"
cat("Saving figure to", filename, "\n")
png(filename, width = 7, height = 3, res = 600, units = "in")
p1 / p2  + plot_layout(heights = c(2, 10))
dev.off()

###########################################################
# DMP function
###########################################################

# helper function to identify DMPs
get_DMPs <- function(design, mVals) {
    
    fit <- lmFit(mVals, design)
    fit2 <- eBayes(fit)
    dmps <- topTable(fit2, coef = 2, number = Inf, adjust.method = "BH", sort.by = "p")
    nsig <- nrow(dmps[dmps$adj.P.Val < 0.1,])
    print(paste("Number of FDR sig DMPs:", nsig))
    return(dmps)
}

###########################################################
# DMP analysis: Responder vs Nonresponder
###########################################################

# responders vs non-responders, control for age and sex
dmp_1 <- get_DMPs(model.matrix(~ responder + Age + Sex, data = targets), mVals)

# responders vs non-responders, control for age
dmp_2 <- get_DMPs(model.matrix(~ responder + Age, data = targets), mVals)

# responders vs non-responders
dmp_3 <- get_DMPs(model.matrix(~ responder, data = targets), mVals)

###########################################################
# DMP analysis: Progression vs Nonprogression
###########################################################

targets$responder <- ifelse(
    targets$Response_RANO %in% c("CR", "PR", "SD"),
    "Responder",
    "NonResponder"
)
targets$responder <- factor(targets$responder)

# responders vs non-responders, control for age and sex
dmp_1 <- get_DMPs(model.matrix(~ responder + Age + Sex, data = targets), mVals)

# responders vs non-responders, control for age
dmp_2 <- get_DMPs(model.matrix(~ responder + Age, data = targets), mVals)

# responders vs non-responders
dmp_3 <- get_DMPs(model.matrix(~ responder, data = targets), mVals)

###########################################################
# DMP analysis: Progression vs Nonprogression
###########################################################

common <- intersect(targets[complete.cases(targets),]$Subject, colnames(mVals))
mVals_sub <- mVals[, common]
targets_sub <- targets[match(common, targets$Subject), ]

# CDKN2, control for age and sex (n=12)
dmp_4 <- get_DMPs(model.matrix(~ CDKN2 + Age + Sex, data = targets_sub), mVals_sub)

# CDKN2, control for age (n=36)
dmp_5 <- get_DMPs(model.matrix(~ CDKN2 + Age, data = targets_sub), mVals_sub)

# CDKN2 (n=95)
dmp_6 <- get_DMPs(model.matrix(~ CDKN2, data = targets_sub), mVals_sub)

###########################################################
# DMR analysis
###########################################################

# get cpg annotation
myannotation <- cpg.annotate(
    datatype = "array",
    object = as.matrix(mVals_sub),
    what = "M",
    analysis.type = "differential",
    design = model.matrix(~ CDKN2, data = targets_sub),
    contrasts = FALSE,
    coef = "responderResponder",
    arraytype = "EPICv1",
    fdr = 0.1
)

# get dmrs
dmrs <- extractRanges(
    dmrcate(myannotation, lambda = 1000, C = 2), genome = "hg19"
)
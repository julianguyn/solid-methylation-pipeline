# trying differential methylation with different versions of minfi normalization

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
bVals <- fread("data/results/data/reprocessing/bVals.tsv", data.table = FALSE)
rownames(bVals) <- bVals$V1
bVals$V1 <- NULL

# get m values
mVals <- fread("data/results/data/reprocessing/mVals.tsv", data.table = FALSE)
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

targets$responder <- factor(meta$response, levels = c("NR", "R"))
targets$Sex <- factor(targets$Sex)

###########################################################
# Target evaluation (not used anymore)
###########################################################

#toPlot <- targets[order(targets$PFS_Days, decreasing = TRUE),]
#toPlot$Subject <- factor(toPlot$Subject, levels = toPlot$Subject)
#toPlot$Response_RANO <- factor(toPlot$Response_RANO, levels = names(mrecist_pal))

#p1 <- ggplot(toPlot, aes(x = Subject, y = 1, fill = CDKN2)) +
#    geom_tile(color = "black") +
#    scale_fill_manual(values = c("#2F4858", "#86BBD8"), na.value = "white") +
#    theme_void()

#p2 <- ggplot(toPlot, aes(x = Subject, y = PFS_Days, fill = Response_RANO)) +
#    geom_bar(stat = "identity") +
#    scale_fill_manual("RECIST", values = mrecist_pal) +
#    theme_minimal() +
#    theme(
#        panel.border = element_rect(),
#        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#        axis.title.x = element_blank()
#    ) +
#    labs(y = "PFS (Days)")

#filename <- "data/results/figures/targets/targets.png"
#cat("Saving figure to", filename, "\n")
#png(filename, width = 7, height = 3, res = 600, units = "in")
#p1 / p2  + plot_layout(heights = c(2, 10))
#dev.off()

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

# responders vs non-responders, control for age and sex (0)
dmp_1 <- get_DMPs(model.matrix(~ responder + Age + Sex, data = targets), mVals)

# responders vs non-responders, control for age
dmp_2 <- get_DMPs(model.matrix(~ responder + Age, data = targets), mVals)

# responders vs non-responders
dmp_3 <- get_DMPs(model.matrix(~ responder, data = targets), mVals)

###########################################################
# DMR analysis
###########################################################

# get cpg annotation
myannotation <- cpg.annotate(
    datatype = "array",
    object = as.matrix(mVals),
    what = "M",
    analysis.type = "differential",
    design = model.matrix(~ responder, data = targets),
    contrasts = FALSE,
    coef = "responderR",
    arraytype = "EPICv1",
    fdr = 0.1
)

# get dmrs
dmrs <- extractRanges(
    dmrcate(myannotation, lambda = 1000, C = 2), genome = "hg19"
)

# https://www.bioconductor.org/packages//release/bioc/vignettes/minfi/inst/doc/minfi.html
# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html

# load libraries
suppressPackageStartupMessages({
    library(stringr)
    library(minfi)
    library(IlluminaHumanMethylationEPICmanifest)
    library(IlluminaHumanMethylationEPICv2manifest)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    library(ggplot2)
    library(ggpubr)
})

###########################################################
# Load in data
###########################################################

# load in sample metadata
meta <- read.csv("metadata/meta_data_research.csv")

# TODO: FIND MAPPING OF METADATA IDS AND IDAT IDS

# load in CpG annotation for both arrays
annEPICv1 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

###########################################################
# Make IDAT meta
###########################################################

# get filenames
targets <- data.frame(
    Sample_Name = sub("_.*", "", list.files("data/rawdata/idats")),
    Basename = sub("_...\\.idat", "", list.files("data/rawdata/idats", full.names = TRUE))
) |> unique()

# get slide and array
targets$Sentrix_ID <- targets$Sentrix_Position <- targets$Array <- NULL
for (i in 1:nrow(targets)) {
    vars <- str_split_1(list.files("data/rawdata/idats")[i], "_")
    targets$Sentrix_ID[i] <- vars[2]
    targets$Sentrix_Position[i] <- vars[3]
}

# annotate array type
targets$Array <- ifelse(
    targets$Sample_Name %in% c("GSM9325970", "GSM9325971", "GSM9325972"),
    "EPIC",     # IlluminaHumanMethylationEPIC
    "EPICv2"    # IlluminaHumanMethylationEPICv2
)

###########################################################
# Get RGSets
###########################################################

# IlluminaHumanMethylationEPIC (3 samples)
rgSet_EPICv1 <- read.metharray.exp(targets = targets[targets$Array == "EPIC",])

# IlluminaHumanMethylationEPICv2 (26 samples)
rgSet_EPICv2 <- read.metharray.exp(targets = targets[targets$Array == "EPICv2",])


###########################################################
# Compare normalization methods
###########################################################

# helper function to get normalized values
compare_norm <- function(RGSet, label) {

    MSet.raw <- preprocessRaw(RGSet)
    MSet.illumina <- preprocessIllumina(RGSet)
    MSet.swan <- preprocessSWAN(RGSet)
    #GRSet.quantile <- preprocessQuantile(RGSet)
    #GRSet.funnorm <- preprocessFunnorm(RGSet)
    GRSet.noob <- preprocessNoob(RGSet)

    p1 <- densityPlot(MSet.raw, main = "Raw", legend = FALSE)
    p2 <- densityPlot(MSet.illumina, main = "Illumina", legend = FALSE)
    p3 <- densityPlot(MSet.swan, main = "SWAN", legend = FALSE)
    #p4 <- densityPlot(GRSet.quantile, main = "Quantile", legend = FALSE)
    #p5 <- densityPlot(GRSet.funnorm, main = "Funnorm", legend = FALSE)
    p6 <- densityPlot(GRSet.noob, main = "Noob", legend = FALSE)

    p7 <- plotQC(getQC(MSet.raw))
    p8 <- plotQC(getQC(MSet.illumina))
    p9 <- plotQC(getQC(MSet.swan))
    #p10 <- plotQC(getQC(GRSet.quantile))
    #p11 <- plotQC(getQC(GRSet.funnorm))
    p12 <- plotQC(getQC(GRSet.noob))

    filename <- paste0("data/results/figures/normalization/", label, ".png")
    png(filename, width = 8, height = 12, res = 600, units = "in")
    print({
        par(mfrow = c(4, 2))
        densityPlot(MSet.raw, main = "Raw", legend = FALSE)
        plotQC(getQC(MSet.raw))
        densityPlot(MSet.illumina, main = "Illumina", legend = FALSE)
        plotQC(getQC(MSet.illumina))
        densityPlot(MSet.swan, main = "SWAN", legend = FALSE)
        plotQC(getQC(MSet.swan))
        densityPlot(GRSet.noob, main = "Noob", legend = FALSE)
        plotQC(getQC(GRSet.noob))
    })
    dev.off()
}

compare_norm(rgSet_EPICv1, label = "EPICv1")
compare_norm(rgSet_EPICv2, label = "EPICv2")

###########################################################
# Probe and Sample QC
###########################################################

# get detection pvals
detPv1 <- detectionP(rgSet_EPICv1)
detPv2 <- detectionP(rgSet_EPICv2)

# plot mean p-values per sample
toPlot <- c(colMeans(detPv1), colMeans(detPv2))
pal <- ifelse(
    toPlot > 0.01,
    "#7C0B2B",
    ifelse(
        names(toPlot) %in% names(colMeans(detPv1)),
        "#545E56", "#917C78")
)
names(toPlot) <- gsub("_.*", "", names(toPlot))
names(pal) <- gsub("_.*", "", names(pal))

filename <- "data/results/figures/qc/detection_pvals.png"
png(filename, width = 6, height = 4, res = 600, units = "in")
barplot(toPlot, las=2, cex.names=0.8, 
        ylab="Mean detection p-values\n",
        col=pal)
abline(h=0.01, col="#7C0B2B", lty=2)
legend("topright", 
       legend = c("EPICv1", "EPICv2", "Failed (detP > 0.01)"),
       fill = c("#545E56", "#917C78", "#7C0B2B"),
       bty = "n")
dev.off()

# set QC thresholds
detP_threshold <- 0.01
sample_failure_rate <- 0.05
probe_failure_rate <- 0.1

# get failed samples
failed_samples_v1 <- colMeans(detPv1 > detP_threshold, na.rm = TRUE) > sample_failure_rate
failed_samples_v2 <- colMeans(detPv2 > detP_threshold, na.rm = TRUE) > sample_failure_rate

# remove failed samples
rgSet_EPICv1 <- rgSet_EPICv1[,!failed_samples_v1]
rgSet_EPICv2 <- rgSet_EPICv2[,!failed_samples_v2] # removing GSM9325977

###########################################################
# Normalization
###########################################################

# normalize with preprocessIllumina()
MSetv1 <- preprocessIllumina(rgSet_EPICv1)
MSetv2 <- preprocessIllumina(rgSet_EPICv2)

# convert to GenomicRatioSet (remove if change to Quantile normalization)
GRSetv1 <- ratioConvert(mapToGenome(MSetv1))
GRSetv2 <- ratioConvert(mapToGenome(MSetv2))

###########################################################
# Remove failed probes
###########################################################

# helper function to remove failed probes
remove_failed_probes <- function(detP, GRSet) {

    # get failed probes
    failed_probes <- rowMeans(detPv1 > detP_threshold, na.rm = TRUE) > probe_failure_rate
    common_probes <- intersect(rownames(GRSet), names(failed_probes))

    # remove failed probes
    to_remove <- failed_probes[common_probes]
    if (length(to_remove) > 0) GRSet <- GRSet[!to_remove]
    cat("Removing", sum(to_remove), "failed probes\n")
    cat("Remaining probes:", nrow(GRSet), "\n")

    return(GRSet)
}

GRSetv1 <- remove_failed_probes(detPv1, GRSetv1)
GRSetv2 <- remove_failed_probes(detPv2, GRSetv2)

###########################################################
# Remove bad probes
###########################################################

# helper function to remove cross reactive and sex chr probes
remove_bad_probes <- function(ann, GRSet) {
    to_remove <- c(
        # cross reactive probes
        rownames(ann)[
            !is.na(ann$Probe_rs)
            # can optionally keep these?
            | !is.na(ann$CpG_rs)
            | !is.na(ann$SBE_rs)
        ],
        # probes on sex chrs
        rownames(ann)[ann$chr %in% c("chrX", "chrY")]
    )
    
    cat("Removing", length(to_remove), "cross reactive or sex chr probes\n")
    if (length(to_remove) > 0) GRSet <- GRSet[!rownames(GRSet) %in% to_remove,]

    # remove probes with SNPs (can optionally drop by maf)
    before <- nrow(GRSet)
    GRSet <- dropLociWithSnps(GRSet)
    cat("Removing", nrow(GRSet) - before, "probes with SNPs\n")
    cat("Remaining probes:", nrow(GRSet), "\n")

    return(GRSet)

}

GRSetv1 <- remove_bad_probes(annEPICv1, GRSetv1)
GRSetv2 <- remove_bad_probes(annEPICv2, GRSetv2)

###########################################################
# Get beta & M values
###########################################################

# combine arrays
GRSet <- combineArrays(GRSet_v1, GRSet_v2, outType = "IlluminaHumanMethylationEPIC")

# get beta values
betas <- getBeta(GRSet)

# get M values
mVals <- getM(GRSet)

###########################################################
# Keep complete probes
###########################################################

# keep probes complete across samples
complete_probes <- complete.cases(mVals) & apply(mVals, 1, function(x) all(is.finite(x)))
cat("Removing", nrow(mVals) - length(complete_probes), "probes\n")

betas <- betas[complete_probes,]
mVals <- mVals[complete_probes,]
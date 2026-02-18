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
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
})


# ---------------------------------------------------------
# Make directories
# ---------------------------------------------------------

outdirs <- c(
    "data/results/figures/normalization",
    "data/results/figures/qc",
    "data/results/data"
)

for (dir in outdirs) {
    if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    }
}

###########################################################
# Load in data
###########################################################

# load in idat and solid id mapping
idat_meta <- read.table("metadata/idat_meta.tsv", header = TRUE)
idat_meta$SolidID[idat_meta$SolidID == "SOLID-006"] <- "SOLID-009"

# load in sample metadata
meta <- read.csv("metadata/meta_data_research.csv")
meta$idat <- idat_meta$SampleID[match(meta$Subject, idat_meta$SolidID)]

# load in CpG annotation for both arrays
annEPICv1 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

###########################################################
# Make targets file
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

# add pheno data
meta <- meta[match(targets$Sample_Name, meta$idat),]
targets <- cbind(targets, meta[,c(1, 3:12, 14:25)])

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

# helper funciton to make density plot
plot_density_beta <- function(GRSet, main, convert = TRUE) {

    if (convert == TRUE) GRSet <- ratioConvert(mapToGenome(GRSet))
    bVals <- getBeta(GRSet)
    test <- reshape2::melt(bVals)

    p <- ggplot(test, aes(x = value, color = Var2)) +
        geom_density() +
        theme_minimal() + 
        theme(
            panel.border = element_rect(linewidth = 1.5),
            legend.position = "none",
            axis.title = element_text(hjust = 0.5)
        ) +
        labs(x = "Beta", y = "Density", title = main)
    return(p)
}

# helper function to get normalized values
compare_norm <- function(RGSet, label) {

    MSet.raw <- preprocessRaw(RGSet)
    MSet.illumina <- preprocessIllumina(RGSet)
    MSet.swan <- preprocessSWAN(RGSet)
    GRSet.quantile <- preprocessQuantile(RGSet)
    GRSet.funnorm <- preprocessFunnorm(RGSet, sex = pData(RGSet)$Sex)
    MSet.noob <- preprocessNoob(RGSet)

    p1 <- plot_density_beta(MSet.raw, main = "Raw", convert = TRUE)
    p2 <- plot_density_beta(MSet.illumina, main = "Illumina", convert = TRUE)
    p3 <- plot_density_beta(MSet.swan, main = "SWAN", convert = TRUE)
    p4 <- plot_density_beta(GRSet.quantile, main = "Quantile", convert = FALSE)
    p5 <- plot_density_beta(GRSet.funnorm, main = "Funnorm", convert = FALSE)
    p6 <- plot_density_beta(MSet.noob, main = "Noob", convert = TRUE)

    filename <- paste0("data/results/figures/normalization/", label, ".png")
    png(filename, width = 12, height = 6, res = 600, units = "in")
    print({ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)})
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
        "#545E56", "#AA9995")
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
       fill = c("#545E56", "#AA9995", "#7C0B2B"),
       bty = "n")
dev.off()

# set QC thresholds
detP_threshold <- 0.01
sample_failure_rate <- 0.05
probe_failure_rate <- 0.1

# get failed samples
failed_samples_v1 <- colMeans(detPv1 > detP_threshold, na.rm = TRUE) > sample_failure_rate
failed_samples_v2 <- colMeans(detPv2 > detP_threshold, na.rm = TRUE) > sample_failure_rate

cat("Number of failed samples:", sum(failed_samples_v1, failed_samples_v2), "\n")
cat("Failed sample:", names(failed_samples_v2[which(failed_samples_v2 == TRUE)]), "\n")

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

    cat("Starting with", nrow(GRSet), "probes\n")

    # get failed probes
    failed_probes <- rowMeans(detPv1 > detP_threshold, na.rm = TRUE) > probe_failure_rate
    common_probes <- intersect(rownames(GRSet), names(failed_probes))

    # remove failed probes
    to_remove <- failed_probes[common_probes]
    if (length(to_remove) > 0) GRSet <- GRSet[!to_remove]
    cat("---Removing", sum(to_remove), "failed probes\n")
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
# Combine arrays
###########################################################

# https://bioconductor.org/packages//release/data/annotation/vignettes/IlluminaHumanMethylationEPICv2manifest/inst/doc/IlluminaHumanMethylationEPICv2manifest.html

# get probes from manifests
probe1 <- getManifestInfo(IlluminaHumanMethylationEPICmanifest, "locusNames")
probe2 <- getManifestInfo(IlluminaHumanMethylationEPICv2manifest, "locusNames")

probe1 <- unique(probe1)
probe2 <- gsub("_.*$", "", probe2)  # remove suffix
probe2 <- unique(probe2)

# identify common probes
common_probes_unfiltered <- intersect(probe1, probe2) # n = 722642
cat(length(common_probes), "common between EPIC and EPICv2\n")

# identify common probes across GRSets
common_probes <- intersect(
    common_probes_unfiltered, 
    intersect(rownames(GRSetv1), gsub("_.*$", "", rownames(GRSetv2)))
) # n = 482264
cat(length(common_probes), "common between GRSets\n")

###########################################################
# Get beta & M values
###########################################################

# combine arrays
# GRSet <- combineArrays(GRSetv1, GRSetv2, outType = "IlluminaHumanMethylationEPIC")

# get beta values
bVals1 <- getBeta(GRSetv1)
bVals2 <- getBeta(GRSetv2)
bVals2 <- do.call(
    rbind,
    tapply(
        1:nrow(bVals2),
        gsub("_.*$", "", rownames(bVals2)),
        function(ind) { colMeans(bVals2[ind, , drop = FALSE]) },
        simplify = FALSE
    ))

# get M values
mVals1 <- getM(GRSetv1)
mVals2 <- getM(GRSetv2)
mVals2 <- do.call(
    rbind,
    tapply(
        1:nrow(mVals2),
        gsub("_.*$", "", rownames(mVals2)),
        function(ind) { colMeans(mVals2[ind, , drop = FALSE]) },
        simplify = FALSE
    ))

# keep common probes
bVals1 <- bVals1[common_probes,]
bVals2 <- bVals2[common_probes,]
mVals1 <- mVals1[common_probes,]
mVals2 <- mVals2[common_probes,]

# merge matrices
bVals <- cbind(bVals1, bVals2) |> as.data.frame()
mVals <- cbind(mVals1, mVals2) |> as.data.frame()

###########################################################
# Keep complete probes
###########################################################

# keep probes complete across samples
complete_probes <- complete.cases(bVals) & apply(bVals, 1, function(x) all(is.finite(x)))
cat("Removing", nrow(bVals) - length(complete_probes), "probes\n")

bVals <- bVals[complete_probes,]
mVals <- mVals[complete_probes,]

###########################################################
# Save matrices
###########################################################

write.table(bVals, file = "data/results/data/bVals.tsv", quote = FALSE, row.names = TRUE)
write.table(mVals, file = "data/results/data/mVals.tsv", quote = FALSE, row.names = TRUE)
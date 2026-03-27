# reprocessing based on manuscript

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

source("workflow/scripts/utils/palettes.R")
source("workflow/scripts/utils/assess_normalization.R")

set.seed(101)

# ---------------------------------------------------------
# Make directories
# ---------------------------------------------------------

outdirs <- c(
    "data/results/data/reprocessing"
)

for (dir in outdirs) {
    if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    }
}

###########################################################
# Load in data
###########################################################

# load in RGSets and targets
load("data/procdata/RGSets_targets.RData")
targets$match <- gsub(".*/", "", targets$Basename)
colnames(rgSet_EPICv1) <- targets$Subject[match(colnames(rgSet_EPICv1), targets$match)]
colnames(rgSet_EPICv2) <- targets$Subject[match(colnames(rgSet_EPICv2), targets$match)]


# load in CpG annotation for both arrays
annEPICv1 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

###########################################################
# Probe and Sample QC
###########################################################

# get detection pvals
detPv1 <- detectionP(rgSet_EPICv1)
detPv2 <- detectionP(rgSet_EPICv2)

# set QC thresholds
detP_threshold <- 0.01
sample_failure_rate <- 0.05
probe_failure_rate <- 0.01

# get failed samples
failed_samples_v1 <- colMeans(detPv1 > detP_threshold, na.rm = TRUE) > sample_failure_rate
failed_samples_v2 <- colMeans(detPv2 > detP_threshold, na.rm = TRUE) > sample_failure_rate

cat("Number of failed samples:", sum(failed_samples_v1, failed_samples_v2), "\n")
cat("Failed sample:", names(failed_samples_v2[which(failed_samples_v2 == TRUE)]), "\n")

# remove failed samples
rgSet_EPICv1 <- rgSet_EPICv1[,!failed_samples_v1]
rgSet_EPICv2 <- rgSet_EPICv2[,!failed_samples_v2] # removing SOLID-040

###########################################################
# Normalization
###########################################################

# normalize with one function only
MSetv1 <- preprocessIllumina(rgSet_EPICv1)
MSetv2 <- preprocessIllumina(rgSet_EPICv2)

# convert to GenomicRatioSet
GRSetv1 <- ratioConvert(mapToGenome(MSetv1))
GRSetv2 <- ratioConvert(mapToGenome(MSetv2))

###########################################################
# Remove failed probes
###########################################################

# helper function to remove failed probes
remove_failed_probes <- function(detP, GRSet, ann) {

    cat("Starting with", nrow(GRSet), "probes\n")

    # get failed probes
    failed_probes <- rowMeans(detPv1 > detP_threshold, na.rm = TRUE) > probe_failure_rate
    common_probes <- intersect(rownames(GRSet), names(failed_probes))

    # remove failed probes
    to_remove <- failed_probes[common_probes]
    if (length(to_remove) > 0) GRSet <- GRSet[!to_remove]
    cat("---Removing", sum(to_remove), "failed probes\n")
    cat("Remaining probes:", nrow(GRSet), "\n")

    # remove bad probes
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

GRSetv1 <- remove_failed_probes(detPv1, GRSetv1, annEPICv1)
GRSetv2 <- remove_failed_probes(detPv2, GRSetv2, annEPICv2)


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
cat(length(common_probes_unfiltered), "common between EPIC and EPICv2\n")

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

write.table(bVals, file = "data/results/data/reprocessing/bVals.tsv", quote = FALSE, row.names = TRUE)
write.table(mVals, file = "data/results/data/reprocessing/mVals.tsv", quote = FALSE, row.names = TRUE)

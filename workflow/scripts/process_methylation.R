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

# load in CpG annotation
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

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
# Get RGsets
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
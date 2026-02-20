# load libraries
suppressPackageStartupMessages({
    library(stringr)
    library(minfi)
    library(IlluminaHumanMethylationEPICmanifest)
    library(IlluminaHumanMethylationEPICv2manifest)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
})

source("workflow/scripts/utils/palettes.R")

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
# Save checkpoint
###########################################################

save(
    targets, rgSet_EPICv1, rgSet_EPICv2,
    file = "data/procdata/RGSets_targets.RData"
)

###########################################################
# Get GRSets
###########################################################

# get GRSets for EPICv1
RGSet <- rgSet_EPICv1
GRSet.raw <- ratioConvert(mapToGenome(preprocessRaw(RGSet)))
GRSet.ill <- ratioConvert(mapToGenome(preprocessIllumina(RGSet)))
GRSet.swn <- ratioConvert(mapToGenome(preprocessSWAN(RGSet)))
GRSet.qnt <- preprocessQuantile(RGSet)
GRSet.fun <- preprocessFunnorm(RGSet, sex = pData(RGSet)$Sex)
GRSet.nob <- ratioConvert(mapToGenome(preprocessNoob(RGSet)))

save(
    GRSet.raw, GRSet.ill, GRSet.swn,
    GRSet.qnt, GRSet.fun, GRSet.nob,
    file = "data/procdata/EPICv1_GRSets.RData"
)

# get GRSets for EPICv2
RGSet <- rgSet_EPICv2
GRSet.raw <- ratioConvert(mapToGenome(preprocessRaw(RGSet)))
GRSet.ill <- ratioConvert(mapToGenome(preprocessIllumina(RGSet)))
GRSet.swn <- ratioConvert(mapToGenome(preprocessSWAN(RGSet)))
GRSet.qnt <- preprocessQuantile(RGSet)
GRSet.fun <- preprocessFunnorm(RGSet, sex = pData(RGSet)$Sex)
GRSet.nob <- ratioConvert(mapToGenome(preprocessNoob(RGSet)))

save(
    GRSet.raw, GRSet.ill, GRSet.swn,
    GRSet.qnt, GRSet.fun, GRSet.nob,
    file = "data/procdata/EPICv2_GRSets.RData"
)
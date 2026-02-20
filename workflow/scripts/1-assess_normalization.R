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
    library(umap)
})

source("workflow/scripts/utils/palettes.R")
source("workflow/scripts/utils/assess_normalization.R")

# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]

valid <- c("EPICv1", "EPICv2")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}

###########################################################
# Load in data
###########################################################

# load in targets
load("data/procdata/RGSets_targets.RData")

# load in GRSets
load(paste0("data/procdata/", analysis, "_GRSets.RData"))

###########################################################
# Plot density plots
###########################################################

p1 <- plot_density_beta(GRSet.raw, main = "Raw")
p2 <- plot_density_beta(GRSet.ill, main = "Illumina")
p3 <- plot_density_beta(GRSet.swn, main = "SWAN")
p4 <- plot_density_beta(GRSet.qnt, main = "Quantile")
p5 <- plot_density_beta(GRSet.fun, main = "Funnorm")
p6 <- plot_density_beta(GRSet.nob, main = "Noob")

filename <- paste0("data/results/figures/normalization/", analysis, "_density.png")
png(filename, width = 12, height = 6, res = 600, units = "in")
print({ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)})
dev.off()

###########################################################
# Plot individual violin plots
###########################################################

# only ran on EPICv2, kept old plots for v1

plot_violin_beta(GRSet.raw, analysis, main = "Raw")
plot_violin_beta(GRSet.ill, analysis, main = "Illumina")
plot_violin_beta(GRSet.swn, analysis, main = "SWAN")
plot_violin_beta(GRSet.qnt, analysis, main = "Quantile")
plot_violin_beta(GRSet.fun, analysis, main = "Funnorm")
plot_violin_beta(GRSet.nob, analysis, main = "Noob")

###########################################################
# Plot umaps
###########################################################

p1 <- plot_umap_mval(GRSet.raw, main = "Raw")
p2 <- plot_umap_mval(GRSet.ill, main = "Illumina")
p3 <- plot_umap_mval(GRSet.swn, main = "SWAN")
p4 <- plot_umap_mval(GRSet.qnt, main = "Quantile")
p5 <- plot_umap_mval(GRSet.fun, main = "Funnorm")
p6 <- plot_umap_mval(GRSet.nob, main = "Noob")

filename <- paste0("data/results/figures/normalization/", analysis, "_umap.png")
png(filename, width = 8, height = 6, res = 600, units = "in")
print({ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)})
dev.off()
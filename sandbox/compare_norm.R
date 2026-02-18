#' Old version of compare_norm()
#' previously used in process_methylation.R

# helper function to get normalized values
compare_norm <- function(RGSet, label) {

    MSet.raw <- preprocessRaw(RGSet)
    MSet.illumina <- preprocessIllumina(RGSet)
    MSet.swan <- preprocessSWAN(RGSet)
    GRSet.quantile <- preprocessQuantile(RGSet)
    GRSet.funnorm <- preprocessFunnorm(RGSet, sex = pData(RGSet)$Sex)
    MSet.noob <- preprocessNoob(RGSet)

    p1 <- densityPlot(MSet.raw, main = "Raw", legend = FALSE)
    p2 <- densityPlot(MSet.illumina, main = "Illumina", legend = FALSE)
    p3 <- densityPlot(MSet.swan, main = "SWAN", legend = FALSE)
    p4 <- densityPlot(GRSet.quantile, main = "Quantile", legend = FALSE)
    p5 <- densityPlot(GRSet.funnorm, main = "Funnorm", legend = FALSE)
    p6 <- densityPlot(MSet.noob, main = "Noob", legend = FALSE)

    p7 <- plotQC(getQC(MSet.raw))
    p8 <- plotQC(getQC(MSet.illumina))
    p9 <- plotQC(getQC(MSet.swan))
    p10 <- plotQC(getQC(GRSet.quantile))
    p11 <- plotQC(getQC(GRSet.funnorm))
    p12 <- plotQC(getQC(MSet.noob))

    filename <- paste0("data/results/figures/normalization/", label, ".png")
    png(filename, width = 12, height = 6, res = 600, units = "in")
    print({
        par(mfrow = c(2, 4))
        densityPlot(MSet.raw, main = "Raw", legend = FALSE)
        densityPlot(MSet.illumina, main = "Illumina", legend = FALSE)
        densityPlot(MSet.swan, main = "SWAN", legend = FALSE)
        densityPlot(GRSet.noob, main = "Noob", legend = FALSE)
        plotQC(getQC(MSet.raw))
        plotQC(getQC(MSet.illumina))
        plotQC(getQC(MSet.swan))
        plotQC(getQC(GRSet.noob))
    })
    dev.off()
}
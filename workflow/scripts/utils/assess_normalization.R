plot_density_beta <- function(GRSet, main) {

    bVals <- getBeta(GRSet)
    toPlot <- reshape2::melt(bVals)
    toPlot$RECIST <- targets$Response_RANO[match(toPlot$Var2, gsub(".*/", "", targets$Basename))]

    p <- ggplot(toPlot, aes(x = value, group = Var2, color = RECIST)) +
        geom_density() +
        scale_color_manual(values = mrecist_pal) +
        theme_minimal() + 
        theme(
            panel.border = element_rect(linewidth = 1.5),
            legend.position = "none",
            axis.title = element_text(hjust = 0.5)
        ) +
        labs(x = "Beta", y = "Density", title = main)
    return(p)
}

plot_violin_beta <- function(GRSet, analysis, main) {

    bVals <- getBeta(GRSet)
    toPlot <- reshape2::melt(bVals)
    toPlot$RECIST <- targets$Response_RANO[match(toPlot$Var2, gsub(".*/", "", targets$Basename))]
    toPlot$ID <- targets$Subject[match(toPlot$Var2, gsub(".*/", "", targets$Basename))]

    p <- ggplot(toPlot, aes(x = ID, y = value, fill = RECIST)) +
        geom_violin() +
        scale_fill_manual(values = mrecist_pal) +
        theme_minimal() + 
        theme(
            panel.border = element_rect(linewidth = 1.5),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none",
            axis.title = element_text(hjust = 0.5)
        ) +
        labs(y = "Beta", title = main)

    width <- switch(
        analysis,
        EPICv1 = 5,
        EPICv2 = 12
    )

    filename <- paste0("data/results/figures/normalization/", analysis, "_", main, "_violin.png")
    png(filename, width = width, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

plot_umap_mval <- function(GRSet, main) {

    mVals <- getM(GRSet)
    toPlot <- reshape2::melt(mVals)

    pheno <- data.frame(
        SampleID = colnames(mVals),
        RECIST = targets$Response_RANO[match(colnames(mVals), gsub(".*/", "", targets$Basename))],
        ID = targets$Subject[match(colnames(mVals), gsub(".*/", "", targets$Basename))]
    )

    # umap
    umap_df <- umap(mVals)$layout
    colnames(umap_df) <- paste0("UMAP", 1:2)
    umap_df <- umap_df[match(pheno$SampleID, rownames(umap_df)),]
    umap_df <- cbind(umap_df, pheno) |> as.data.frame()

    p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = RECIST)) +
        geom_point(shape = 21, size = 2) +
        scale_fill_manual(values = mrecist_pal) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            legend.key.size = unit(0.5, 'cm')
        ) +
        labs(title = main)
    return(p)
}
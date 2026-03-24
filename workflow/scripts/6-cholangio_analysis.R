# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(rGREAT)
    library(ggplot2)
    library(viridis)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

dmrs <- fread("data/rawdata/cholangio/cholangio_baseline_IDHwt_Vs_IDHmut_DiffExpr_limma_Result_p_0.1.txt", data.table = FALSE)
dmrs <- dmrs[!is.na(dmrs$logFC),] # 2933 dmrs
hyper <- dmrs[which(dmrs$logFC > 0),] # 1302 dmrs
hypo <- dmrs[which(dmrs$logFC < 0),] # 1631 dmrs

###########################################################
# Create genomic ranges
###########################################################

# get coordinates
hyper_coords <- do.call(rbind, strsplit(hyper$Region_ID, "\\.")) |> as.data.frame()
hypo_coords <- do.call(rbind, strsplit(hypo$Region_ID, "\\.")) |> as.data.frame()

# make genomic ranges object
hyper_gr <- GRanges(seqnames = hyper_coords$V1, ranges = IRanges(as.numeric(hyper_coords$V2), as.numeric(hyper_coords$V3)))
hypo_gr <- GRanges(seqnames = hypo_coords$V1, ranges = IRanges(as.numeric(hypo_coords$V2), as.numeric(hypo_coords$V3)))

###########################################################
# GREAT analysis
###########################################################

run_great <- function(gr, label) {
    job <- submitGreatJob(gr, species = "hg19", genome = "hg19", help = FALSE)
    tbl <- getEnrichmentTables(job)

    # save results
    mf <- as.data.frame(tbl[1])
    bp <- as.data.frame(tbl[2])
    cc <- as.data.frame(tbl[3])

    # save column names
    cols <- gsub("GO.Molecular.Function.", "", colnames(mf))
    colnames(mf) <- colnames(bp) <- colnames(cc) <- cols

    mf$Label <- "MF"
    bp$Label <- "BP"
    cc$Label <- "CC"

    # combine results and filter
    res <- rbind(bp, mf, cc)
    res <- res[res$Hyper_Adjp_BH < 0.05,]
    res <- res[order(res$Hyper_Adjp_BH),]

    filename <- paste0("data/results/data/cholangio/", label, ".tsv")
    write.table(res, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
    return(res)
}

hyper_res <- run_great(hyper_gr, "hypermethylation_GREAT") # 244 pathways
hypo_res <- run_great(hypo_gr, "hypomethylation_GREAT") # 172

###########################################################
# Plot GREAT analysis
###########################################################

plot_GREAT <- function(great, label) {

    great <- great[order(great$Hyper_Fold_Enrichment, decreasing = TRUE),]
    toPlot <- great[1:30,]
    toPlot$name <- factor(toPlot$name, levels=rev(toPlot$name))

    p <- ggplot(toPlot, aes(x = name, y = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH)) +
        geom_point(aes(size = Hyper_Total_Genes, shape = Label)) +
        geom_segment(aes(x = name, xend = name, y = 0, yend = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH), size = 1) +
        coord_cartesian(clip = "off") + 
        coord_flip() +
        scale_alpha(range = c(1, 0.2)) +
        guides(
            shape = guide_legend(title = "Ontology", override.aes = list(size = 4)),
            color = guide_legend(override.aes = list(shape = 19, size = 4))
        ) +
        scale_color_viridis_c(option = "rocket", direction = -1, end = 0.85) +
        theme_classic() + 
        theme(
            legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
            legend.position = c(0.95, 0.05),   # bottom-right
            legend.justification = c("right", "bottom"),
            legend.background = element_rect(fill = alpha('white', 0.7), color = NA),
            plot.title = element_text(hjust = 0.5, size = 12), 
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        ) + 
        labs(x = "GO Term", y = "Fold Enrichment", size = "Genes", color = "Adjusted\nP-Value", title = label)

    filename <- paste0("data/results/figures/cholangio/", label, "_GREAT.png")
    png(filename, width = 8, height = 7, res = 600, units = "in")
    print(p)
    dev.off()
}

plot_GREAT(hyper_res, "hypermethylation")
plot_GREAT(hypo_res, "hypomethylation")

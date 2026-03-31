# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(rGREAT)
    library(ggplot2)
    library(viridis)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(pheatmap)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

# load in dmrs
dmrs <- fread("data/rawdata/cholangio/cholangio_baseline_IDHwt_Vs_IDHmut_DiffExpr_limma_Result_p_0.1.txt", data.table = FALSE)
dmrs <- dmrs[!is.na(dmrs$logFC),] # 2933 dmrs
hyper <- dmrs[which(dmrs$logFC > 0),] # 1302 dmrs
hypo <- dmrs[which(dmrs$logFC < 0),] # 1631 dmrs

# load in matrix
mat <- readRDS("data/rawdata/cholangio/MeDIP_Solid_cholangiocarcinoma_log2CPM.rds")

# load in metadata
pheno <- fread("data/rawdata/cholangio/Pheno_cholangio_all_IDHwt_Vs_IDHmut.txt")

# load in annotation
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

###########################################################
# Quick heatmap
###########################################################

toPlot <- mat[rownames(mat) %in% dmrs$Region_ID,]

p <- pheatmap(
    toPlot,
    scale = "row",
    cluster_rows = TRUE,
    cutree_rows = 2,
    treeheight_row = 20,
    show_rownames = FALSE
)

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
# Annotate genomic regions
###########################################################

# annotate regions
hyper_anno <- annotatePeak(hyper_gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")
hypo_anno <- annotatePeak(hypo_gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")

# create toPlot
hyper <- hyper@annoStat
hypo <- hypo@annoStat
hyper$Label <- "Hyper"
hypo$Label <- "Hypo"
toPlot <- rbind(hyper, hypo)

# plot gene features
filename <- paste0("data/results/figures/cholangio/peakanno.png")
png(filename, width = 5, height = 5, res = 600, units = "in")
ggplot(toPlot, aes(fill=Feature, y=Frequency, x=Label)) + 
    geom_bar(position="fill", stat="identity", color = "black") +
    scale_fill_manual(values = genefeat_pal) +
    theme_minimal() + 
    labs(y = "Percentage (%)")
dev.off()

###########################################################
# Keep promoter regions
###########################################################

# 133 hypermethylated promoters
hyper_df <- as.data.frame(hyper_anno)
hyper_promoter <- hyper_df[grep("Promoter", hyper_df$annotation), ]

# 225 hypomethylated promoters
hypo_df <- as.data.frame(hypo_anno)
hypo_promoter <- hypo_df[grep("Promoter", hypo_df$annotation), ]

# make genomic ranges object
hyper_gr <- GRanges(seqnames = hyper_promoter$seqnames, ranges = IRanges(hyper_promoter$start, hyper_promoter$end))
hypo_gr <- GRanges(seqnames = hypo_promoter$seqnames, ranges = IRanges(hypo_promoter$start, hypo_promoter$end))


###########################################################
# Run GO analysis
###########################################################

# run GO
run_go <- function(df, label) {

    genes <- unique(df$geneId)
    ontologies <- c("BP", "MF", "CC")

    # run enrichGO
    go_list <- lapply(ontologies, function(o) {
        enrichGO(
            gene = genes,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = o,
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05,
            readable = TRUE
        ) |> as.data.frame()
    })
    names(go_list) <- ontologies

    # bind results
    bp <- go_list$BP
    mf <- go_list$MF
    cc <- go_list$CC

    if (nrow(bp) > 0) bp$Label <- "BP"
    if (nrow(mf) > 0) mf$Label <- "MF"
    if (nrow(cc) > 0) cc$Label <- "CC"

    res <- rbind(bp, mf, cc)
    filename <- paste0("data/results/data/cholangio/", label, "_GO.tsv")
    write.table(res, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
    return(res)
}

# genes from all regions
hyper_go_all <- run_go(hyper_df, "hyper_all") # 11 results
hypo_go_all <- run_go(hypo_df, "hypo_all") # 64 results

# genes from promoter regions only
hyper_go_promoter <- run_go(hyper_promoter, "hyper_promoter") # no results
hypo_go_promoter <- run_go(hypo_promoter, "hypo_promoter") # 7 results


###########################################################
# Plot GO analysis
###########################################################

plot_go <- function(go_res, label) {

    go_res <- go_res[order(go_res$p.adjust, decreasing = TRUE),]
    go_res$Description <- factor(go_res$Description, levels = go_res$Description)

    if (label == "hypo_promoter") {
        w <- 5
        h <- 3.5
    } else {
        w <- 7
        h <- 5
        go_res <- go_res[1:10,]
    }

    # add gene ratio
    go_res$GeneRatio <- sapply(go_res$GeneRatio, function(x) {
    vals <- strsplit(x, "/")[[1]]
    as.numeric(vals[1]) / as.numeric(vals[2])
    })

    filename <- paste0("data/results/figures/cholangio/", label, "_GO.png")
    png(filename, width = 5, height = 3.5, res = 600, units = "in")
    print(ggplot(go_res, aes(x = Label, y = Description, color = p.adjust, size = GeneRatio)) +
        geom_point() +
        scale_color_viridis_c(option = "rocket", direction = -1, begin = 0.15, end = 0.85) +
        theme_bw() +
        theme(axis.title.y = element_blank()) +
        labs(x = "Ontology"))
    dev.off()

}

plot_go(hyper_go_all, "hyper_all")
plot_go(hypo_go_all, "hypo_all")
plot_go(hypo_go_promoter, "hypo_promoter")

###########################################################
# GREAT analysis
###########################################################

# helper function to run GREAT anlaysis
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

hyper_res <- run_great(hyper_gr, "hyper_promoter_GREAT") # 244 all, 0 promoter
hypo_res <- run_great(hypo_gr, "hypo_promoter_GREAT") # 172 all, 0 promoter

###########################################################
# Plot GREAT analysis
###########################################################

# helper function to plot GREAT analysis results
plot_GREAT <- function(great, label) {

    great <- great[great$Hyper_Total_Genes > 10,]
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
        scale_color_viridis_c(option = "rocket", direction = -1, begin = 0.15, end = 0.85) +
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
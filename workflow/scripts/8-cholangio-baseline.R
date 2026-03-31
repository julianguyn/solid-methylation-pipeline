# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(limma)
    library(ggplot2)
    library(ggpubr)
    library(RColorBrewer)
    library(viridis)
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(rGREAT)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

mat <- readRDS("data/rawdata/cholangio/MeDIP_Solid_cholangiocarcinoma_log2CPM.rds")

pheno <- fread("data/rawdata/cholangio/Pheno_cholangio_all_IDHwt_Vs_IDHmut.txt") #64
#pheno2 <- fread("data/rawdata/cholangio/Pheno_cholangio_baseline_IDHwt_Vs_IDHmut.txt") #16

# keep only cholangio samples
mat <- mat[,colnames(mat) %in% pheno$MeDIP_ID]
pheno <- pheno[match(colnames(mat), pheno$MeDIP_ID),] # 38

ids <- unique(pheno$solid_patient)
for (id in ids) {
    subset <- pheno$response[pheno$solid_patient == id]
    cat(id, ":", unique(subset), "\n")
}

# load in annotation
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

###########################################################
# Get variables
###########################################################

# rank date for each sample
pheno$date <- as.Date(pheno$blood_collection_date, format = "%Y %b %d")
pheno <- pheno[order(pheno$solid_patient, pheno$date), ]
pheno$timepoint <- unlist(lapply(as.numeric(table(pheno$solid_patient)), seq_len))
pheno$timepoint <- factor(paste0("T", pheno$timepoint), levels = paste0("T", 1:12))

# add final response
baseline <- pheno[pheno$response == "baseline",]
baseline$final_response <- ifelse(
    baseline$solid_patient %in% pheno$solid_patient[pheno$response == "SD"],
    "SD", "PD"
)

###########################################################
# Keep baseline samples
###########################################################

baseline <- baseline[baseline$phase == "c1d1",]
baseline <- unique(baseline[,c(2, 6, 10, 15, 18)])
# pheno[,c(2,6, 10, 17)]
# baseline[,c(2, 6, 10, 17, 18)]
# keep only c1d1 to be consistent (remove screens)

# subset matrix
mat <- mat[,match(baseline$MeDIP_ID, colnames(mat))]

###########################################################
# Run linear model (SD-PD)
###########################################################

# set variables
labels <- factor(baseline$final_response)

# model
design <- model.matrix(~ labels) # 465 DMRs with FDR<0.05
fit <- lmFit(mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.05)

write.csv(res, file = "data/results/data/cholangio/baseline_sd_vs_pd.csv", quote = FALSE, row.names = TRUE)

#ggplot() +
#    geom_point(data = res, aes(x = logFC, y = -log10(adj.P.Val)), color = "gray") +
#    geom_point(data = res[res$logFC > 2 & res$adj.P.Val < 0.05,], aes(x = logFC, y = -log10(adj.P.Val)), color = "red") +
#    geom_point(data = res[res$logFC < -2 & res$adj.P.Val < 0.05,], aes(x = logFC, y = -log10(adj.P.Val)), color = "blue") +
#    theme_bw()

hyper <- res[res$logFC > 2 & res$adj.P.Val < 0.05,] #54
hypo <- res[res$logFC < -2 & res$adj.P.Val < 0.05,] #411

###########################################################
# Create genomic ranges
###########################################################

# get coordinates
hyper_coords <- do.call(rbind, strsplit(rownames(hyper), "\\.")) |> as.data.frame()
hypo_coords <- do.call(rbind, strsplit(rownames(hypo), "\\.")) |> as.data.frame()

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
hyper <- hyper_anno@annoStat
hypo <- hypo_anno@annoStat
hyper$Label <- "Hyper"
hypo$Label <- "Hypo"
toPlot <- rbind(hyper, hypo)
toPlot$Feature <- factor(toPlot$Feature, levels = names(genefeat_pal))

# plot gene features
filename <- paste0("data/results/figures/cholangio/baseline_sdpd_peakanno.png")
png(filename, width = 4, height = 4, res = 600, units = "in")
ggplot(toPlot, aes(fill=Feature, y=Frequency, x=Label)) + 
    geom_bar(position="fill", stat="identity", color = "black") +
    scale_fill_manual(values = genefeat_pal) +
    theme_minimal() + 
    labs(y = "Percentage (%)")
dev.off()

###########################################################
# Keep promoter regions
###########################################################

# 11 hypermethylated promoters
hyper_df <- as.data.frame(hyper_anno)
hyper_promoter <- hyper_df[grep("Promoter", hyper_df$annotation), ]

# 38 hypomethylated promoters
hypo_df <- as.data.frame(hypo_anno)
hypo_promoter <- hypo_df[grep("Promoter", hypo_df$annotation), ]

# make genomic ranges object
hyper_gr_promoter <- GRanges(seqnames = hyper_promoter$seqnames, ranges = IRanges(hyper_promoter$start, hyper_promoter$end))
hypo_gr_promoter <- GRanges(seqnames = hypo_promoter$seqnames, ranges = IRanges(hypo_promoter$start, hypo_promoter$end))

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
            pvalueCutoff = 0.1,
            qvalueCutoff = 0.1,
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
    filename <- paste0("data/results/data/cholangio/baseline_sdpd_", label, "_GO.tsv")
    write.table(res, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
    print(nrow(res))
    return(res)
}

# genes from all regions
hyper_go_all <- run_go(hyper_df, "hyper_all") # 0 results
hypo_go_all <- run_go(hypo_df, "hypo_all") # 1 result

# genes from promoter regions only
hyper_go_promoter <- run_go(hyper_promoter, "hyper_promoter") # 99 results
hypo_go_promoter <- run_go(hypo_promoter, "hypo_promoter") # 11 results

###########################################################
# Plot GO analysis
###########################################################

plot_go <- function(go_res, label) {

    go_res <- go_res[order(go_res$p.adjust, decreasing = TRUE),]
    go_res$Description <- factor(go_res$Description, levels = go_res$Description)

    if (label == "hypo_promoter") {
        w <- 6
        h <- 3.5
    } else {
        w <- 6.5
        h <- 4
        go_res <- go_res[1:15,]
    }

    # add gene ratio
    go_res$GeneRatio <- sapply(go_res$GeneRatio, function(x) {
    vals <- strsplit(x, "/")[[1]]
    as.numeric(vals[1]) / as.numeric(vals[2])
    })

    filename <- paste0("data/results/figures/cholangio/baseline_sdpd_", label, "_GO.png")
    png(filename, width = w, height = h, res = 600, units = "in")
    print(ggplot(go_res, aes(x = Label, y = Description, color = p.adjust, size = GeneRatio)) +
        geom_point() +
        scale_color_viridis_c(option = "rocket", direction = -1, begin = 0.15, end = 0.85) +
        theme_bw() +
        theme(axis.title.y = element_blank(), legend.key.size = unit(0.5, 'cm')) +
        labs(x = "Ontology"))
    dev.off()

}

plot_go(hyper_go_promoter, "hyper_promoter")
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
    res <- res[res$Hyper_Adjp_BH < 0.1,]
    res <- res[order(res$Hyper_Adjp_BH),]

    filename <- paste0("data/results/data/cholangio/baseline_sdpd_", label, ".tsv")
    write.table(res, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
    print(nrow(res))
    return(res)
}

hyper_res <- run_great(hyper_gr, "hyper_GREAT") # 0
hyper_res_promoter <- run_great(hyper_gr_promoter, "hyper_promoter_GREAT") #0
hypo_res <- run_great(hypo_gr, "hypo_GREAT") # 7 all
hypo_res_promoter <- run_great(hypo_gr_promoter, "hypo_promoter_GREAT") #0

###########################################################
# Plot GREAT analysis
###########################################################

# helper function to plot GREAT analysis results
plot_GREAT <- function(great, label) {

    #great <- great[great$Hyper_Total_Genes > 10,]
    toPlot <- great[order(great$Hyper_Fold_Enrichment, decreasing = TRUE),]
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
            legend.key.size = unit(0.3, 'cm'), text = element_text(size = 10),
            #legend.position = c(0.95, 0.05),   # bottom-right
            #legend.justification = c("right", "bottom"),
            #legend.background = element_rect(fill = alpha('white', 0.7), color = NA),
            plot.title = element_text(hjust = 0.5, size = 12), 
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        ) + 
        labs(x = "GO Term", y = "Fold Enrichment", size = "Genes", color = "Adjusted\nP-Value", title = label)

    filename <- paste0("data/results/figures/cholangio/baseline_sdpd_", label, "_GREAT.png")
    png(filename, width = 10, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

plot_GREAT(hypo_res, "hypomethylation")
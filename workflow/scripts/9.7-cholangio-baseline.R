# try pooling all samples based on baseline labels

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
    library(ComplexHeatmap)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

mat <- readRDS("data/rawdata/cholangio/MeDIP_Solid_cholangiocarcinoma_log2CPM.rds")
pheno <- fread("data/rawdata/cholangio/Pheno_cholangio_all_IDHwt_Vs_IDHmut.txt") #64

# keep only cholangio samples
common <- intersect(pheno$MeDIP_ID, colnames(mat))
mat <- mat[,match(common, colnames(mat))]
pheno <- pheno[match(common, pheno$MeDIP_ID),] #56 samples from 10 patients

# load in annotation
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

###########################################################
# Get variables
###########################################################

# add final response
pheno$final_response <- ifelse(
    pheno$solid_patient %in% pheno$solid_patient[pheno$response == "SD"],
    "SD", "PD"
)
# recode Solid-10
pheno$final_response[pheno$solid_patient == "Solid-10"] <- "PD"
unique(pheno[,c(2, 16)][order(pheno$final_response),])

# remove PD sample from SD patients
sd_patients <- pheno$solid_patient[pheno$final_response == "SD"]
pheno <- pheno[-which(pheno$solid_patient %in% sd_patients & pheno$response == "PD"),]

# keep common samples
common <- intersect(pheno$MeDIP_ID, colnames(mat))
mat <- mat[,match(common, colnames(mat))]
pheno <- pheno[match(common, pheno$MeDIP_ID),] #51 samples from 10 patients


###########################################################
# Run linear model (SD-PD)
###########################################################

# set variables
labels <- factor(pheno$final_response)
patientid <- factor(pheno$solid_patient)

# model
design <- model.matrix(~ labels + patientid) # 1140 DMRs with FDR<0.05
fit <- lmFit(mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.05)
 
write.csv(res, file = "data/results/data/cholangio/frombaseline_all_sd_vs_pd.csv", quote = FALSE, row.names = TRUE)


hyper <- res[res$logFC > 0 & res$adj.P.Val < 0.05,] #279
hypo <- res[res$logFC < 0 & res$adj.P.Val < 0.05,] #861

###########################################################
# Create genomic ranges
###########################################################

# get coordinates
hyper_coords <- do.call(rbind, strsplit(rownames(hyper), "\\.")) |> as.data.frame()
hypo_coords <- do.call(rbind, strsplit(rownames(hypo), "\\.")) |> as.data.frame()
bg_coords <- do.call(rbind, strsplit(rownames(res), "\\.")) |> as.data.frame()

# make genomic ranges object
hyper_gr <- GRanges(seqnames = hyper_coords$V1, ranges = IRanges(as.numeric(hyper_coords$V2), as.numeric(hyper_coords$V3)))
hypo_gr <- GRanges(seqnames = hypo_coords$V1, ranges = IRanges(as.numeric(hypo_coords$V2), as.numeric(hypo_coords$V3)))
bg_gr <- GRanges(seqnames = bg_coords$V1, ranges = IRanges(as.numeric(bg_coords$V2), as.numeric(bg_coords$V3)))

###########################################################
# Annotate genomic regions
###########################################################

# annotate regions
hyper_anno <- annotatePeak(hyper_gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")
hypo_anno <- annotatePeak(hypo_gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")
bg_anno <- annotatePeak(bg_gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")

# create toPlot
hyper <- hyper_anno@annoStat
hypo <- hypo_anno@annoStat
bg <- bg_anno@annoStat
hyper$Label <- "Hyper"
hypo$Label <- "Hypo"
bg$Label <- "Background"
toPlot <- rbind(hyper, hypo, bg)
toPlot$Feature <- factor(toPlot$Feature, levels = names(genefeat_pal))

# plot gene features
filename <- paste0("data/results/figures/cholangio/FINALrmPD_all_sdpd_peakanno.png")
png(filename, width = 5, height = 4, res = 600, units = "in")
ggplot(toPlot, aes(fill=Feature, y=Frequency, x=Label)) + 
    geom_bar(position="fill", stat="identity", color = "black") +
    scale_fill_manual(values = genefeat_pal) +
    theme_minimal() + 
    labs(y = "Percentage (%)")
dev.off()

###########################################################
# Keep promoter regions
###########################################################

# 42 hypermethylated promoters
hyper_df <- as.data.frame(hyper_anno)
hyper_promoter <- hyper_df[grep("Promoter", hyper_df$annotation), ]

# 122 hypomethylated promoters
hypo_df <- as.data.frame(hypo_anno)
hypo_promoter <- hypo_df[grep("Promoter", hypo_df$annotation), ]

bg_df <- as.data.frame(bg_anno)
bg_promoter <- bg_df[grep("Promoter", bg_df$annotation),]

# make genomic ranges object
hyper_gr_promoter <- GRanges(seqnames = hyper_promoter$seqnames, ranges = IRanges(hyper_promoter$start, hyper_promoter$end))
hypo_gr_promoter <- GRanges(seqnames = hypo_promoter$seqnames, ranges = IRanges(hypo_promoter$start, hypo_promoter$end))
bg_gr_promoter <- GRanges(seqnames = bg_promoter$seqnames, ranges = IRanges(bg_promoter$start, bg_promoter$end))

###########################################################
# Run GO analysis
###########################################################

# run GO
run_go <- function(df, bg_gr, label) {

    genes <- unique(df$geneId)
    ontologies <- c("BP", "MF", "CC")

    # run enrichGO
    go_list <- lapply(ontologies, function(o) {
        enrichGO(
            gene = genes,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = o,
            universe = unique(bg_gr$geneID),
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
    filename <- paste0("data/results/data/cholangio/frombaseline_all_sdpd_", label, "_GO.tsv")
    write.table(res, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
    print(nrow(res))
    return(res)
}

# genes from all regions
hyper_go_all <- run_go(hyper_df, bg_gr, "hyper_all") # 0 results
hypo_go_all <- run_go(hypo_df, bg_gr, "hypo_all") # 6 result

# genes from promoter regions only
hyper_go_promoter <- run_go(hyper_promoter, bg_gr_promoter, "hyper_promoter") # 8 results
hypo_go_promoter <- run_go(hypo_promoter, bg_gr_promoter, "hypo_promoter") # 33 results

###########################################################
# Plot GO analysis
###########################################################

plot_go <- function(go_res, label) {

    go_res <- go_res[order(go_res$p.adjust, decreasing = TRUE),]
    go_res$Description <- factor(go_res$Description, levels = go_res$Description)

    if (label == "hypo_promoter") {
        w <- 5
        h <- 3.5
        go_res <- go_res[1:15,]
        print(go_res$Description)
    } else {
        w <- 4.5
        h <- 3.5
    }

    # add gene ratio
    go_res$GeneRatio <- sapply(go_res$GeneRatio, function(x) {
    vals <- strsplit(x, "/")[[1]]
    as.numeric(vals[1]) / as.numeric(vals[2])
    })

    filename <- paste0("data/results/figures/cholangio/FINAL_all_sdpd_", label, "_GO.png")
    png(filename, width = w, height = h, res = 600, units = "in")
    print(ggplot(go_res, aes(x = Label, y = Description, color = p.adjust, size = GeneRatio)) +
        geom_point() +
        scale_color_viridis_c(option = "rocket", direction = -1, begin = 0.15, end = 0.85) +
        theme_bw() +
        theme(axis.title.y = element_blank(), legend.key.size = unit(0.5, 'cm')) +
        guides(
            size = guide_legend(order = 1),
            color = guide_colorbar(order = 2)
        ) +
        labs(x = "Ontology", color = "FDR"))
    dev.off()

}

plot_go(hypo_go_promoter, "hypo_promoter")

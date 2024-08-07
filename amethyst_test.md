```bash
singularity shell ~/singularity/amethyst/
```

```R
library(amethyst)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(plyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)

#download files
system("mkdir Downloads")
download.file("https://adeylabopen.s3.us-west-2.amazonaws.com/amethyst/pbmc_vignette.h5", "./Downloads/pbmc_vignette.h5", method = "curl") # Contains site-level methylation information for each cell
download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/pbmc_vignette/pbmc_vignette_cellInfo.txt", "./Downloads/pbmc_vignette_cellInfo.txt") # Summary QC statistics for each cell
download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/pbmc_vignette/pbmc_vignette.annot", "./Downloads/pbmc_vignette.annot") # Simulated batch metadata

obj <- createObject()
obj <- addCellInfo(obj, file = "./Downloads/pbmc_vignette_cellInfo.txt")
obj <- addAnnot(obj, file = "./Downloads/pbmc_vignette.annot", name = "batch") 
#annot is cellID<tab>metadatacolumn
head(obj@metadata)
plt<-ggplot(obj@metadata, aes(x = cov)) + geom_histogram(bins = 10) 
ggsave(file="./Downloads/cov_plot.pdf",plt)

#filter cells
obj@metadata <- obj@metadata |> dplyr::filter(cov > 100000 & cov < 40000000)

#set up h5 paths (basically same format as an annot)
obj@h5paths <- data.frame(row.names = rownames(obj@metadata), paths = rep("./Downloads/pbmc_vignette.h5", length(rownames(obj@metadata))))
head(obj@h5paths)

obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = 1) 

#run on 100kb windows
obj@genomeMatrices[["cg_100k_score"]] <- makeWindows(obj,
                                                     stepsize = 100000, 
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 1, 
                                                     index = "chr_cg", 
                                                     nmin = 2, 
                                                     species = "human") 

#filter windows by cell coverage
obj@genomeMatrices[["cg_100k_score"]] <- obj@genomeMatrices[["cg_100k_score"]][rowSums(!is.na(obj@genomeMatrices[["cg_100k_score"]])) >= 45, ] #45 is because there are 50 cells

dimEstimate(obj, genomeMatrices = c("cg_100k_score"), dims = c(10), threshold = 0.95)
#cg_100k_score 
#            7 
set.seed(111)
obj@reductions[["irlba"]] <- runIrlba(obj, genomeMatrices = c("cg_100k_score"), dims = c(7), replaceNA = c(0))
obj@reductions[["irlba_regressed"]] <- regressCovBias(obj, reduction = "irlba") # Optional; helps reduce coverage 
obj <- runCluster(obj, k_phenograph = 10, reduction = "irlba_regressed") # consider increasing k_phenograph to 50 for larger datasets
obj <- runUmap(obj, neighbors = 5, dist = 0.05, method = "euclidean", reduction = "irlba_regressed") 
obj <- runTsne(obj, perplexity = 10, method = "euclidean", theta = 0.2, reduction = "irlba_regressed") 

p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle("UMAP")
p2 <- dimFeature(obj, colorBy = cluster_id, reduction = "tsne") + ggtitle("TSNE")
plt<-plot_grid(p1, p2)
ggsave(plt,file="./Downloads/umap.pdf")

#run on windows (bed format)
download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/pbmc_vignette/pbmc_highconfidence_dmrs.bed", "./Downloads/pbmc_dmr.bed") 
obj@genomeMatrices[["pbmc_dmrs"]] <- makeWindows(obj, bed = "./Downloads/pbmc_dmr.bed", type = "CG", metric = "percent", threads = 1, index = "chr_cg", nmin = 2, species = "human")

obj@genomeMatrices[["pbmc_dmrs"]] <- obj@genomeMatrices[["pbmc_dmrs"]][rowSums(!is.na(obj@genomeMatrices[["pbmc_dmrs"]])) >= 10, ] 
dimEstimate(obj, genomeMatrices = c("pbmc_dmrs"), dims = c(10), threshold = 0.90)
set.seed(111)
obj@reductions[["irlba"]] <- runIrlba(obj, genomeMatrices = c("pbmc_dmrs"), dims = c(8), replaceNA = c(0))
obj@reductions[["irlba_regressed"]] <- regressCovBias(obj, reduction = "irlba")
obj <- runCluster(obj, k_phenograph = 10, reduction = "irlba_regressed") # consider increasing k_phenograph to 50 for larger datasets
obj <- runUmap(obj, neighbors = 5, dist = 0.05, method = "euclidean", reduction = "irlba_regressed") 
plt1<-dimFeature(obj, colorBy = cluster_id, reduction = "umap")
plt2<-dimFeature(obj, colorBy = cluster_id) + facet_wrap(vars(batch)) # Batch is simulated to illustrate function utility. Any column in the metadata will work.
plt<-plot_grid(plt1, plt2)
ggsave(plt,file="./Downloads/umap_dmr.pdf")

#plot by feature
plt<-sampleComp(obj, groupBy = "batch", colorBy = "cluster_id") 
ggsave(plt,file="./Downloads/sample_comp.pdf")

#plot by coverage and %mCG
p1 <- dimFeature(obj, colorBy = log(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
p2 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
plt<-plot_grid(p1, p2)
ggsave(plt,file="./Downloads/sample_cov_umap.pdf")


#annotate data
obj@ref <- makeRef("hg38")
cluster1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 1,
                                         step = 1000,
                                         smooth = 3,
                                         species = "human",
                                         index = "chr_cg",
                                         groupBy = "cluster_id",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj@genomeMatrices[["cg_cluster_tracks"]] <- cluster1kbwindows[["pct_matrix"]]

plt<-heatMap(obj, 
        genes = c("SPI1", "CD2", "S100A8",  "CD79A", "CD3G", "ELANE", "MPO", 
           "S100A8", "MPEG1", "IRF8", "CD74", "GZMK", "CD3E", "CD3D", "KIR2DL4", "KLRB1"), 
        matrix = "cg_cluster_tracks", 
        nrow = 4,
        legend = FALSE,
        width = 1000)
ggsave(plt,file="./Downloads/marker_heatmap.pdf")

plt<-histograM(obj, genes = c("ELANE", "MPEG1", "SPI1", "CD2", "CD3D"), matrix = "cg_cluster_tracks", legend = F, width = 1000)
ggsave(plt,file="./Downloads/marker_histogram.pdf")


#plot gene promoters as dotplots
protein_coding <- unique(obj@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding" & seqid != "chrM") |> dplyr::pull(gene_name))
obj@genomeMatrices[["cg_promoter"]] <- makeWindows(obj, 
                                                     genes = protein_coding,
                                                     promoter = TRUE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 1, 
                                                     index = "chr_cg", 
                                                     nmin = 2, 
                                                     species = "human") 
# subsetting to genes with values in at least 5 cells
obj@genomeMatrices[["cg_promoter"]] <- obj@genomeMatrices[["cg_promoter"]][rowSums(!is.na(obj@genomeMatrices[["cg_promoter"]])) >= 10, ]

genes <- c("SPI1", "CD19", "CD2", "CD6", "CD8A", "CD4", "CSF1R", "GATA1", "CD79A", "CD3G", "ELANE", "MPO", "ITGAM",
           "S100A8", "MPEG1", "FN1", "IRF8", "CD74", "RORA", "GZMK", "CD3E", "CD3D", "MEIS1", "KIR2DL4")
plt<-dotM(obj, genes = genes, groupBy = "cluster_id", matrix = "cg_promoter") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "#999999")) + scale_size(range = c(1, 8))
ggsave(plt,file="./Downloads/marker_dotplot.pdf")


#find cluster markers
cluster_promoter_markers <- findClusterMarkers(obj, 
                                               matrix = "cg_promoter", 
                                               genes = genes, 
                                               threads = 1)
cluster_promoter_markers <- cluster_promoter_markers |> dplyr::filter(p.adj < 0.05) # Not many results because it's a small dataset
head(cluster_promoter_markers)
plt<-dotM(obj, genes = cluster_promoter_markers$gene, groupBy = "cluster_id", matrix = "cg_promoter") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "#999999")) + scale_size(range = c(1, 12))
ggsave(plt,file="./Downloads/denovomarkers_dotplot.pdf")


#DMR
pal <- c("#F9AB60", "#E7576E", "#630661", "#B5DCA5") # makePalette(option = 7, n = 4) 

dmrs <- testDMR(cluster1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 5, nminGroup = 5) # or use cluster1kbwindows[["sum_matrix"]] and rename
dmrs <- filterDMR(dmrs, method = "bonferroni", filter = TRUE, pThreshold = 0.01, logThreshold = 2)
head(dmrs)
collapsed_dmrs <- collapseDMR(obj, dmrs, maxDist = 4000, minLength = 2000, reduce = T, annotate = T) 
head(collapsed_dmrs)
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(test, direction) |> dplyr::summarise(n = n()), 
       aes(y = test, x = n, fill = test)) + geom_col() + 
  facet_grid(vars(direction), scales = "free_y") + scale_fill_manual(values = pal) + theme_classic()
ggsave(plt,file="./Downloads/met_per_dmr.pdf")

top_dmrs <- collapsed_dmrs |> 
  dplyr::group_by(type, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(test, direction) |> slice_min(n = 1, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)

heatMap(obj, matrix = "cg_type_tracks", regions = top_dmrs$location[top_dmrs$direction == "hypo"], nrow = 2, legend = FALSE, width = 1000, arrowOverhang = 0)

top_dmrs <- collapsed_dmrs |> 
  dplyr::group_by(test, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(test, direction) |> slice_min(n = 1, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)

plt<-heatMap(obj, matrix = "cg_type_tracks", regions = top_dmrs$location[top_dmrs$direction == "hypo"], nrow = 2, legend = FALSE, width = 1000, arrowOverhang = 0)
ggsave(plt,file="./Downloads/dmr_denovo_markers.pdf")
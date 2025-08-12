library(Seurat)
library(glue)
library(patchwork)
library(ggplot2)
library(purrr)

plot_cell_cycle_distribution <- function(sample_id, study) {

  seu <- readRDS(glue("/dataVolume/storage/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_cnv_seu.rds"))
  DimPlot(seu, group.by = "Phase")

  cell_cycle_phase_per_cluster <- janitor::tabyl(seu@meta.data, Phase, gene_snn_res.0.2) %>%
    tidyr::pivot_longer(-c("Phase"), names_to = "cluster", values_to = "num_cell") %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(percent_cell = num_cell/sum(num_cell))

  # percent cell  ------------------------------
  cell_cycle_distribution = ggplot(cell_cycle_phase_per_cluster, aes(percent_cell, cluster, fill = Phase)) +
    geom_col() +
    theme_bw() +
    labs(title = sample_id, y = "Percent cells")

  phase_dimplot <- DimPlot(seu, group.by = "Phase")

  cluster_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2")

  patchwork <- cell_cycle_distribution / (phase_dimplot + cluster_dimplot)

}

# field ------------------------------
sample_ids <- as.character(glue("SRR1796048{0:4}"))
field_cell_cycle_patchworks <- map(sample_ids, plot_cell_cycle_distribution, "field")

pdf("results/cellcycle_patchworks.pdf")
field_cell_cycle_patchworks
dev.off()

# yang ------------------------------
sample_ids <- as.character(glue("SRR148005{34:43}"))
sample_ids <- sample_ids[c(1:4, 6:10)]
yang_cell_cycle_patchworks <- map(sample_ids, plot_cell_cycle_distribution, "yang")

pdf("../yang_et_al_proj/results/cellcycle_patchworks.pdf")
yang_cell_cycle_patchworks
dev.off()

# wu ------------------------------
sample_ids <- as.character(glue("SRR1388424{0:9}"))
wu_cell_cycle_patchworks <- map(sample_ids, plot_cell_cycle_distribution, "collin")

pdf("../collin_et_al_proj/results/cellcycle_patchworks.pdf")
wu_cell_cycle_patchworks
dev.off()



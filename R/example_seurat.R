
source("R/packages.R")
source("R/functions.R")

# load cellranger dirs
cellranger_dirs <-
  dir_ls("~/myprojects/HughesC/250721_0725Bio-05_HughC_10x3v4_GEX/CellRanger", 
         glob = "*filtered_feature_bc_matrix*",
         recurse = TRUE,
         type = "directory") |> 
    # set_names(str_extract("/CellRanger/([^/]+)/outs/")) |> 
    identity()

# metadata <- 
  dir_ls("~/myprojects/HughesC/250721_0725Bio-05_HughC_10x3v4_GEX", glob = "*meta*", 
         recurse = TRUE)

# test0 <- 
names(cellranger_dirs) <- 
  str_extract(cellranger_dirs, pattern = "(?<=/CellRanger/)[^/]+(?=/outs/)")

sample_mats <- map(cellranger_dirs, Read10X)

seus <- imap(sample_mats, 
            ~Seurat::CreateSeuratObject(
              counts = .x, project = .y, 
            ))

seu1 <- Read10X(cellranger_dirs[[1]])
seu1 <- CreateSeuratObject(seu1, project = names(cellranger_dirs)[[1]])

seu_p <- seurat_process(seu1)

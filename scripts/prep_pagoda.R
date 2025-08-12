
args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
	eval(parse(text = args[[i]]))
}

library(pagoda2)
library(glue)
library(numbat)

countMatrix <- Seurat::Read10X(matrix_dir)

# all cells ------------------------------

## all basic pagoda2 processing with basicP2proc()
pagoda <- basicP2proc(countMatrix, n.cores=as.integer(ncores), min.cells.per.gene=10, 
														n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)

## calculate pathway overdispersion for human
# pagoda <- extendedP2proc(pagoda, organism = 'hs')

saveRDS(pagoda, outrds)

# # create app
# 
# ## create app object
# p2app <- webP2proc(pagoda$p2, title = 'Quick pagoda2 app', go.env = pagoda$go.env)
# 
# # save web app
# saveRDS(p2app, webapp)

# # # aneuploid cells ------------------------------
# 
# nb = Numbat$new(out_dir = numbat_outdir)
# 
# aneuploid_cells <- nb$clone_post %>% 
# 	dplyr::filter(compartment_opt == "tumor") %>% 
# 	dplyr::pull(cell) %>% 
# 	identity()
# 
# aneuploid_countMatrix <- countMatrix[,aneuploid_cells]
# 
# ## all basic pagoda2 processing with basicP2proc()
# aneuploid_pagoda <- basicP2proc(aneuploid_countMatrix, n.cores=ncores, min.cells.per.gene=10,
# 											n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
# 
# ## calculate pathway overdispersion for human
# aneuploid_pagoda <- extendedP2proc(aneuploid_pagoda, organism = 'hs')
# 
# saveRDS(aneuploid_pagoda, aneuploid_outrds)

# create aneuploid app

# ## create app object
# aneuploid_p2app <- webP2proc(aneuploid_pagoda$p2, title = 'Quick pagoda2 app', go.env = aneuploid_pagoda$go.env)
# 
# # save web app
# saveRDS(aneuploid_p2app, aneuploid_webapp)

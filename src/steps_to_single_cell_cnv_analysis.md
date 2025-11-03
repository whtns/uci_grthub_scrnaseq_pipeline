flow chart

priority samples 
  wu: SRR13884245
  field: SRR17960480, SRR17960481, SRR17960482, SRR17960484, 

2. annotate clusters by cell type (cellpypes) [annotate_cellpypes.R]
3. exclude cell types [annotate_cellpypes.R]
  * immune (B2M)
  * rod (NRL)
1. define genotypes from infercnv [manual]
3. record clusters with SCNAs [plotting?]
4. plot cell cycle phase by cluster
5. differential expression between cells w/ and w/o SCNAs not part of cell type outgroups
6. define genotype proportion by cluster

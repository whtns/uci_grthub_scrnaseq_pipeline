#!/usr/bin/env bash

for i in {1..22}
do
  chrom_ending="chr""$i"
  
  cellsnp_vcf = outputdir + "/numbat/{sample}/pileup/{sample}/cellSNP.base.vcf",

  vcfTarget=$(echo "$1" | sed 's/.vcf.gz/_/g')"$chrom_ending".vcf.gz
  
  bcfTarget=$chrom_ending".genotypes.bcf"
  
  geneticMapFile=$2
  
  numThreads=$3
  
done

# eagle --numThreads 6 --vcfTarget output/numbat/SRR13633760/phasing/SRR13633760_chr1.vcf.gz --vcfRef /dataVolume/storage/Homo_sapiens/numbat/1000G_hg38/chr1.genotypes.bcf --geneticMapFile=/usr/local/bin/TOOLS/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz --outPrefix output/numbat/SRR13633760/phasing/SRR13633760_chr1.phased

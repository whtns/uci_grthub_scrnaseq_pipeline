#!/usr/bin/env Rscript 

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

samples = "SRR13633760"
outdir = "output/numbat/SRR13633760/"
label = samples


library(logger, quietly = T)
library(glue, quietly = T)
library(stringr, quietly = T)
library(argparse, quietly = T)
library(data.table, quietly = T)
library(dplyr, quietly = T)
library(vcfR, quietly = T)
library(Matrix, quietly = T)
library(numbat)

samples = str_split(samples, ',')[[1]]
bams = str_split(bams, ',')[[1]]
barcodes = str_split(barcodes, ',')[[1]]
n_samples = length(samples)
genome = ifelse(str_detect(gmap, 'hg19'), 'hg19', 'hg38')
message(paste0('Using genome version: ', genome))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (sample in samples) {
    dir.create(glue('{outdir}/pileup'), showWarnings = FALSE)
    dir.create(glue('{outdir}/phasing'), showWarnings = FALSE)
    dir.create(glue('{outdir}/pileup/{sample}'), showWarnings = FALSE)
}

## pileup

cmds = c()

if (smartseq) {

    cmd = glue(
            'cellsnp-lite', 
            '-S {bams}',
            '-i {barcodes}',
            '-O {outdir}/pileup/{samples}',
            '-R {snpvcf}', 
            '-p {ncores}',
            '--minMAF 0',
            '--minCOUNT 2',
            '--UMItag None',
            '--cellTAG None',
            .sep = ' ')

    cmds = c(cmd)

} else {
    
    for (i in 1:n_samples) {
        
        cmd = glue(
            'cellsnp-lite', 
            '-s {bams[i]}',
            '-b {barcodes[i]}',
            '-O {outdir}/pileup/{samples[i]}',
            '-R {snpvcf}', 
            '-p {ncores}',
            '--minMAF 0',
            '--minCOUNT 2',
            '--UMItag {UMItag}',
            '--cellTAG {cellTAG}',
            .sep = ' ')

        cmds = c(cmds, cmd)

    }

}

cat('Running pileup\n')

script = glue('{outdir}/run_pileup.sh')

list(cmds) %>% fwrite(script, sep = '\n')

system(glue('chmod +x {script}'))

# tryCatch({
#     system(glue('sh {script} &> {outdir}/pileup.log'), intern = TRUE)
# },
# warning = function(w){
#     stop('Pileup failed')
# })
# 
# ## VCF creation
# cat('Creating VCFs\n')
# vcfs = lapply(samples, function(sample){vcfR::read.vcfR(glue('{outdir}/pileup/{sample}/cellSNP.base.vcf'), verbose = F)})
# 
# genotype(label, samples, vcfs, glue('{outdir}/phasing'))
# 
# ## phasing
# eagle_cmd = function(chr, sample) {
#     paste(eagle, 
#         glue('--numThreads {ncores}'), 
#         glue('--vcfTarget {outdir}/phasing/{label}_chr{chr}.vcf.gz'), 
#         glue('--vcfRef {paneldir}/chr{chr}.genotypes.bcf'), 
#         glue('--geneticMapFile={gmap}'), 
#         glue('--outPrefix {outdir}/phasing/{label}_chr{chr}.phased'),
#     sep = ' ')
# }
# 
# cmds = c("#!/usr/bin/env bash")
# 
# for (sample in samples) {
#     cmds = c(cmds, lapply(1:22, function(chr){eagle_cmd(chr, sample)}))
# }
# 
# script = glue('{outdir}/run_phasing.sh')
# 
# list(cmds) %>% fwrite(script, sep = '\n')
# 
# system(glue('chmod 777 {script}'))
# 
# system2(script, stdout = glue("{outdir}/phasing.log"))
# 
# tryCatch({
#     system2(script, stdout = glue("{outdir}/phasing.log"))
#     # system(script)
# },
# warning = function(w){
#     stop('Phasing failed')
# })
# 
# ## Generate allele count dataframe
# cat('Generating allele count dataframes\n')

if (genome == 'hg19') {
    gtf_transcript = gtf_hg19
} else {
    gtf_transcript = gtf_hg38
}

for (sample in samples) {
    print(sample)
    
    # read in phased VCF
    vcf_phased = lapply(1:22, function(chr) {
            vcf_file = glue('{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')
            if (file.exists(vcf_file)) {
                fread(vcf_file) %>%
                    rename(CHROM = `#CHROM`) %>%
                    mutate(CHROM = str_remove(CHROM, 'chr'))   
            }
        }) %>%
        Reduce(rbind, .) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM)))

    pu_dir = glue('{outdir}/pileup/{sample}')

    # pileup VCF
    vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf')) %>% rename(CHROM = `#CHROM`)

    # count matrices
    AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
    DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))

    cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)

    df = preprocess_allele(
        sample = label,
        vcf_pu = vcf_pu,
        vcf_phased = vcf_phased,
        AD = AD,
        DP = DP,
        barcodes = cell_barcodes,
        gtf_transcript = gtf_transcript
    )
    
    fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
    
}

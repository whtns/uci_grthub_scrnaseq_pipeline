
rule pileup_and_phase:
	input:
		bam = outputdir + "cellranger/{sample}/outs/possorted_genome_bam.bam",
		barcodes = "data/{sample}/barcodes.tsv.gz"
		script = "scripts/pileup_and_phase.R"
	output:
		outputdir + "numbat/{sample}/{sample}_allele_counts.tsv.gz"
	log:
		outputdir + "logs/numbat_{sample}.log"
	benchmark:
		outputdir + "benchmarks/numbat_{sample}.txt"
	params:
		gmap = config["gmap"]
		snpvcf = config["snpvcf"]
		paneldir = config["paneldir"]
	conda:
		"../envs/environment_R.yaml"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args label='{wildcards.sample}' samples='{wildcards.sample}' '''
		'''bams='{input.bam}' barcodes='{input.barcodes} gmap='{params.gmap}' snpvcf='{params.snpvcf}' '''
		'''paneldir='{params.paneldir}' outdir='outputdir + "numbat/{wildcards.sample}' ncores='4'" {input.script} {log}'''

# rule numbat:
# 	input:
# 		allele_table = outputdir + "numbat/{sample}/{sample}_allele_counts.tsv.gz",
# 		expression_matrix = count_mat_ATC2,
# 		script = "scripts/run_numbat.R"
# 	output:
# 		outdir = outputdir + "numbat/{sample}/{sample}_allele_counts.tsv.gz"
# 	log:
# 		outputdir + "logs/numbat_{sample}.log"
# 	benchmark:
# 		outputdir + "benchmarks/numbat_{sample}.txt"
# 	params:
# 		gmap = config["gmap"]
# 		snpvcf = config["snpvcf"]
# 		paneldir = config["paneldir"]
# 	conda:
# 		"../envs/environment_R.yaml"
# 	shell:
# 		'''{Rbin} CMD BATCH --no-restore --no-save "--args df_allele_ATC2='{input.allele_table}' count_mat_ATC2='{input.expression_matrix}' label='{wildcards.sample}' ncores='4'" {input.script} {log}'''

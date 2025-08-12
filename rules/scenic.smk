

rule pileup_and_phase:
	input:
		bam_file = "output/STAR/SRR13633760/SRR13633760_Aligned.sortedByCoord.out.bam"
		barcodes = "data/9_4months_Retinoblastoma/barcodes.tsv.gz"
		script = "scripts/pileup_and_phase.R"
	output:
		outdir = outputdir + "numbat/pileup/{sample}"
	log:
		outputdir + "Rout/scenic.Rout"
	benchmark:
		outputdir + "benchmarks/scenic.txt"
	params:
		gmap = config["gmap"]
		snpvcf = config["snpvcf"]
		paneldir = config["paneldir"]
		numbat_dir = outputdir + "numbat/"
	conda:
		"../envs/environment_R.yaml"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args label='test' samples='test' '''
		'''bams='{input.bam_file}' barcodes='{input.barcodes} gmap='{params.gmap}' snpvcf='test' '''
		'''paneldir='{params.paneldir}' outdir='{params.numbat_dir}' ncores='4'" {input.script} {log}'''

rule runscenic:
	input:
		loom_file = outputdir + "scenic/unfiltered.loom"
	output:
		final_loom = outputdir + "scenic/unfiltered-final.loom"
	log:
		outputdir + "Rout/scenic.Rout"
	benchmark:
		outputdir + "benchmarks/scenic.txt"
	params:
		TFs = config['TFs'],
		motifs = config['motifs'],
		feather_db = config['feather_db']
	shell:
	  "nextflow run aertslab/SCENICprotocol -profile docker "
	  "--loom_input {input.loom_file} --loom_output {output.final_loom} "
	  "--TFs {params.TFs} --motifs {params.motifs} --db {params.feather_db} "
	  "--thr_min_genes 1"

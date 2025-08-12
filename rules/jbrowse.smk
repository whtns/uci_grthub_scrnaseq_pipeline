
## ------------------------------------------------------------------------------------ ##
## configure jbrowse
## ------------------------------------------------------------------------------------ ##
## configure jbrowse

rule jbrowsemeta:
  input:
    metatsv = config['metatsv']
  output:
    metatsv = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + config['metatsv']
  threads:
    config['ncores']
  conda:
    "../envs/environment.yaml"
  shell:
    "cp {input.metatsv} {output.metatsv}"

rule jbrowse:
	input:
	  hisatbam = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
	  hisatbai = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
	  hisatbigwig = outputdir + "HISAT2bigwig/{sample}_Aligned.sortedByCoord.out.bw"
	output:
	  bam_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/samples/{sample}.bam",
	  bai_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/samples/{sample}.bam.bai",
	  bigwig_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/samples/{sample}.bw"
	threads:
		config["ncores"]
	params:
	  proj_name = os.path.basename(proj_dir)
	conda:
		"../envs/environment.yaml"
	log:
		outputdir + "logs/jbrowse_{sample}.log"
	shell:
	  "ln -s {input.hisatbam} {output.bam_symlink}; "
	  "ln -s {input.hisatbai} {output.bai_symlink}; "
	  "ln -s {input.hisatbigwig} {output.bigwig_symlink}"

rule jbrowsetracklist:
	input:
		script = "scripts/format_tracklist_json.py",
		refdir = config["refdir"],
		gff = config["gff"],
		gff_tbi = config["gff_tbi"],
		refseq = config["refseq"]
	output:
	  tracklist_json = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/trackList.json",
	  refdir_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/reference/",
	  gff_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/Homo_sapiens.GRCh38.87.sorted.gff3.gz",
	  gff_tbi_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/Homo_sapiens.GRCh38.87.sorted.gff3.gz.tbi",
	  refseq_symlink = "/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/seq/refSeqs.json"
	threads:
		config["ncores"]
	log:
		outputdir + "logs/browsetracklist.log"
	params:
	  proj_name = os.path.basename(proj_dir),
	  metatsv = os.path.basename(proj_dir) + "/" + config['metatsv']
	conda:
		"../envs/environment.yaml"
	shell:
	  "{input.script} '{params.proj_name}' '{params.metatsv}'; "
	  "ln -sr {input.refdir} {output.refdir_symlink}; "
	  "ln -s {input.gff} {output.gff_symlink}; "
	  "ln -s {input.gff_tbi} {output.gff_tbi_symlink}; "
	  "ln -s {input.refseq} {output.refseq_symlink}; "

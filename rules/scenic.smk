
# Rule to subsample loom file to 10% of cells using loompy (Python)
rule subsample_loom:
    input:
        loom_file = f"{OUTPUT_DIR}/loom/{{sample}}.loom"
    output:
        subsampled_loom = f"{OUTPUT_DIR}/loom/{{sample}}_subsampled.loom"
    resources:
        mem_mb = 4000,
        cpus = 2
    log:
        f"{OUTPUT_DIR}/Rout/{{sample}}/subsample_loom.log"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/{{sample}}_subsample_loom.txt"
    threads: 2
    conda:
        "../envs/environment.yaml"
    shell:
        """
        python scripts/subsample_loom.py \
            {input.loom_file} \
            {output.subsampled_loom} \
            --percentage 0.3
        """

# SCENIC rule for running the SCENIC analysis pipeline
rule runscenic:
    input:
        loom_file = f"{OUTPUT_DIR}/loom/{{sample}}_subsampled.loom"
    output:
        scenic_loom = f"{OUTPUT_DIR}/scenic/{{sample}}.loom"
    resources:
        mem_mb = 16000,  # 16GB for SCENIC
        cpus = 8,  
        partition = "standard"
    log:
        f"{OUTPUT_DIR}/Rout/{{sample}}/scenic.Rout"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/{{sample}}_scenic.txt"
    threads: 8
    params:
        TFs = config['TFs'],
        motifs = config['motifs'],
        feather_db = config['feather_db']
    shell:
        """
        module load singularity/3.11.3
        nextflow run aertslab/SCENICprotocol -profile singularity \
        --loom_input {input.loom_file} --loom_output {output.scenic_loom} \
        --TFs {params.TFs} --motifs {params.motifs} --db {params.feather_db} \
        --thr_min_genes 10
        module unload singularity/3.11.3
        """


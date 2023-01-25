rule fastqc:
    input:
        "data/bisulfite-seq/{species}/{sample}.fastq.gz"
    output:
        html="output/qc/fastqc/{species}/{sample}_fastqc.html",
        zip="output/qc/fastqc/{species}/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log: "logs/fastqc/{species}/{sample}.log"
    params:
        outdir = lambda wildcards: "output/qc/fastqc/{}".format(wildcards.species)
    threads: 4
    conda: "../envs/QC.yaml"
    shell: 
        "fastqc --quiet -t {threads} --outdir {params.outdir} {input}  > {log} 2>&1"

rule fastqc_trimmomatic:
    input:
        "output/trimmed/{species}/{sample}_trimmed_{read}.fq.gz"
    output:
        html="output/qc/fastqc/trimmed/{species}/{sample}_trimmed_{read}_fastqc.html",
        zip="output/qc/fastqc/trimmed/{species}/{sample}_trimmed_{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log: "logs/fastqc/{species}/{sample}_trimmed_{read}.log"
    params:
        outdir = lambda wildcards: "output/qc/fastqc/trimmed/{}".format(wildcards.species)
    threads: 4
    conda: "../envs/QC.yaml"
    shell: 
        "fastqc --quiet -t {threads} --outdir {params.outdir} {input}  > {log} 2>&1"

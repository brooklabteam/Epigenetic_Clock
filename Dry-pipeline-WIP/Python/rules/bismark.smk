import os


def getgenome_bismark_genome_prep(wildcards):
    genome = "data/genomes/{s}/{g}/{g}.fa".format(s=wildcards.species, g = config["genome"][wildcards.species])
    #print(genome)
    return genome

rule bismark_genome_prep:
    input: getgenome_bismark_genome_prep
    output: directory("data/genomes/{species}/{genome}/Bisulfite_Genome")
    params:
        genomedir = lambda wildcards: "data/genomes/{species}/{genome}".format(species=wildcards.species, genome = wildcards.genome)
    log: "logs/bismark_genome_prep/{species}/{genome}/bisulfite.log"
    wildcard_constraints:
    	species = "[A-Z]_[a-z]+(Probe)?"
    threads: 40
    conda: "../envs/environment.yml"
    shell:
        """
        bismark_genome_preparation --bowtie2 --verbose {params.genomedir} 2> {log}
        """

rule bismark:
    input:
        read1 = "output/trimmed/{species}/{sample}_trimmed_R1.fq.gz",
        read2 = "output/trimmed/{species}/{sample}_trimmed_R2.fq.gz",
        bisulfiteGenome = "data/genomes/{species}/{genome}/Bisulfite_Genome"
    output: temp("output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.bam"),
    log:
        log1="logs/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.log",
        log2="logs/bismark/{species}/{genome}_{sample}_trimmed_bismark_PE_report.txt",
        log3="logs/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nucleotide_stats.txt"
    params:
        # Params used by Griffin et al 2021
        extra = " -N 1 --non_directional",
        genomedir = lambda wildcards: "data/genomes/{s}/{g}".format(s=wildcards.species, g=wildcards.genome),
        outdir = lambda wildcards: "output/bismark/{s}/".format(s=wildcards.species)
    threads: 40
    conda: "../envs/environment.yml"
    shell:
        """
        ROOTPROJDIR="$(pwd -P)"
        bismark {params.extra} --bowtie2 -p {threads} --nucleotide_coverage {params.genomedir} -1 {input.read1} -2 {input.read2} --basename {wildcards.genome}_{wildcards.sample}_trimmed_bismark --output_dir {params.outdir} 2> {log.log1}
        ln -sf $ROOTPROJDIR/{params.outdir}/{wildcards.genome}_{wildcards.sample}_trimmed_bismark_PE_report.txt {log.log2}
        ln -sf $ROOTPROJDIR/{params.outdir}/{wildcards.genome}_{wildcards.sample}_trimmed_bismark_pe.nucleotide_stats.txt {log.log3}
        """

# bismark function filter_non_conversion (option --threshold 11)

rule bismark_filter_non_conversion:
    input:
        bam="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.bam"
    output:
        bam=temp("output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.bam"),
        other="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_removed_seqs.bam"
    params:
        outdir=lambda wildcards: "output/bismark/{species}/".format(species=wildcards.species),
        threshold = 11 #default 3
    log:
        "logs/bismark/{species}/{genome}_{sample}_filter_non_conversion.log"
    conda: "../envs/environment.yml"
    shell:
        """
        filter_non_conversion --threshold {params.threshold} --paired {input.bam} 2> {log}
        """


rule bismark_deduplicate:
    input:
        bam="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.bam"
    output:
        bam=temp("output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.bam")
    params:
        outdir=lambda wildcards: "output/bismark/{species}/".format(species=wildcards.species)
    log:
        "logs/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.log"
    conda: "../envs/environment.yml"
    shell:
        """
        deduplicate_bismark --paired --bam {input.bam} --output_dir {params.outdir} 2> {log}
        """


rule sort_individual:
    input:
        bam="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.bam"
    output:
        sort="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.bam",
        index="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.bam.bai"
    conda: "../envs/environment.yml"
    shell:
        """
        samtools sort {input.bam} > {output.sort}
        samtools index {output.sort}
        """


rule bamCoverage_individual:
    input:
        bam="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.bam",
        bai="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.bam.bai"
    output:
        "output/bamCoverage/{species}/{genome}_{sample}_bismark.bw"
    params:
        binsize=config['binSize'],
        extra=""
    threads: 8
    conda: "../envs/environment.yml"
    shell:
        """
        bamCoverage {params.extra} -b {input.bam} -o {output}  --binSize {params.binsize} -p {threads}
        """


rule methylation_extractor:
   input:
       bam="output/bismark/{species}/{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.bam"
   output:
       CHH="output/methylation_extracted/{species}/CHH_context_{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.txt.gz",
       CHG="output/methylation_extracted/{species}/CHG_context_{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.txt.gz",
       CpG="output/methylation_extracted/{species}/CpG_context_{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.txt.gz",
   params:
       genomedir = lambda wildcards: "data/genomes/{s}".format(s=wildcards.species),
       outdir = lambda wildcards: "output/methylation_extracted/{s}".format(s=wildcards.species),
       extra = "--comprehensive --report --buffer_size 8G --no_overlap"
   threads: 4
   conda: "../envs/environment.yml"
   log:
       "logs/methylation_extractor/{species}/{genome}_{sample}.log",
   shell:
       """
       bismark_methylation_extractor {params.extra} --gzip --paired-end --multicore {threads}  --genome_folder {params.genomedir} -s {input.bam} --output {params.outdir} 2> {log}
       """

rule bismark2bedGraph:
    input:
        "output/methylation_extracted/{species}/{context}_context_{genome}_{sample}_trimmed_bismark_pe.nonCG_filtered.deduplicated.sorted.txt.gz",
    output:
        bedGraph = "output/bismark_bedgraph/{species}/{genome}_{sample}_{context}.gz",
        cov = "output/bismark_bedgraph/{species}/{genome}_{sample}_{context}.gz.bismark.cov.gz",
    params:
        outdir = lambda wildcards: "output/bismark_bedgraph/{s}".format(s=wildcards.species),
        extra = ""
    conda: "../envs/environment.yml"
    log:
        "logs/bismarkd2bedGraph/{species}/{genome}_{sample}_{context}.log",
    shell:
        """
        bismark2bedGraph {params.extra} --CX {input} -o {wildcards.genome}_{wildcards.sample}_{wildcards.context} --dir {params.outdir} 2> {log}
        """

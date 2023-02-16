import pandas as pd

configfile: "config.yaml"

include: "rules/trimmomatic.smk"
include: "rules/bismark.smk"
include: "rules/fastQC.smk"
include: "rules/QCGraphs.smk"


def get_all_methylation(wildcards):
    samples = pd.read_table("samples.tsv", index_col= False)[["Sample","Species"]]
    CHH_file_path = "output/methylation_extracted/{Species}/CHH_context_{genome}_{Sample}_trimmed_bismark_pe.deduplicated.sorted.txt.gz"
    files = [CHH_file_path.format(genome = r.Species[0], **r) for i,r in samples.iterrows()]
    #print(files)
    if len(files) == 0:
        raise Exception("No samples found! Check your samples file and try again!")
    else:
        return files


rule all_methylation_extracted:
    input: get_all_methylation

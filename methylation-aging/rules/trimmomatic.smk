import pandas as pd


def get_all_trimmomatic_pe(wildcards):
    file_path = "output/trimmed/{Species}/{Sample}_trimmed_{Read}.fq.gz"
    samples = pd.read_table("samples.tsv")[["Species", "Sample", "Read"]]
    input_samples = [file_path.format(**r) for i, r in samples.iterrows()]
    #input_samples = pd.read_table("samples.tsv").Path.to_list()
    print(input_samples)
    if len(input_samples) == 0:
        raise Exception("No samples found! Check pepsamples.tsv and try again!")
    else:
        return input_samples
       

rule all_trimmomatic_pe:
    input: get_all_trimmomatic_pe



def get_trimmomatic_pe(wildcards):
    samples = pd.read_table("samples.tsv", index_col= "Sample")
    sample_rec = samples.loc[wildcards.sample]
    input_samples = {"r1": sample_rec[sample_rec["Read"] == "R1"].Path[0], "r2": sample_rec[sample_rec["Read"] == "R2"].Path[0]}
    #print(input_samples)
    if len(input_samples) == 0:
        raise Exception("No samples found for species {}. Check pepsamples.tsv and try again!".format(wildcards.species))
    else:
        return input_samples

rule trimmomatic_pe:
    input: unpack(get_trimmomatic_pe)
    output:
        r1="output/trimmed/{species}/{sample}_trimmed_R1.fq.gz",
        r2="output/trimmed/{species}/{sample}_trimmed_R2.fq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/{species}/{sample}_trimmed_R1_unpaired.fq.gz",
        r2_unpaired="output/trimmed/{species}/{sample}_trimmed_R2_unpaired.fq.gz"
    log:
        "logs/trimmomatic/{species}/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:data/trimmomatic-adapters/TruSeq3-PE-2.fa:2:40:15", 
                 "SLIDINGWINDOW:5:20"
                ],
        # optional parameters
        extra="",
        compression_level="-9"
    threads: 4
    wrapper:
        "0.74.0/bio/trimmomatic/pe"


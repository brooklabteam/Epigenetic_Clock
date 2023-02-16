configfile: "config.yaml"

def get_all_qc(wildcards):
    import pandas as pd
    import os
    samples = pd.read_table("samples.tsv", index_col= False) #[["Filename","Species"]]
    samples.Filename = [i.rsplit(".")[0] for i in samples.Filename]
    qc_file_path = "output/qc/fastqc/{Species}/{Filename}_fastqc.html"
    qc_trimmed_file_path = "output/qc/fastqc/trimmed/{Species}/{Sample}_trimmed_{Read}_fastqc.html"
    files = [qc_file_path.format(**r) for i,r in samples[["Filename", "Species"]].iterrows()]
    files += [qc_trimmed_file_path.format(**r) for i,r in samples[["Species","Sample","Read"]].iterrows()]
    if len(files) == 0:
        raise Exception("No samples found! Check your samples file and try again!")
    else:
        return files
        
        
rule multiqc:
    input: get_all_qc
    output: "output/qc/multiqc/multiqc.html"
    params:
        input_dir = "output/qc/fastqc",
        output_dir = "output/qc/multiqc"
    conda: "../envs/QC.yaml"
    shell: "multiqc --force --interactive -o {params.output_dir} -n multiqc.html --profile-runtime  {params.input_dir}"

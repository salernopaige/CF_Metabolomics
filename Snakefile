# All the important snakemake options are recorded in cluster_params/config.yaml
# Usage:
# snakemake --profile slurm

import re
from pathlib import Path
import logging

configfile: "config.yaml"

indir = Path(config["indir"])
outdir = Path(config["outdir"])
if "rare_indir" in config:
    rare_indir=Path(config["rare_indir"])
else:
    rare_indir=None

if config["samples_sub"] is not None:
    SAMPLES = config["samples_sub"]
elif config["samples_dir_regex"] is not None:
    SAMPLES = [x.name for x in Path(indir).glob(config["samples_dir_regex"]) if x.is_dir()]
else:
    SAMPLES = [x.name for x in Path(indir).glob("*") if x.is_dir()]

if config["samples_ignore"] is not None:
    SAMPLES = [x for x in SAMPLES if x not in config["samples_ignore"]]

assert len(SAMPLES) != 0, "Could not find any samples"

rule all:
    input:
        expand("%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir), sample=SAMPLES),

def input_reads(wildcards):
    if config["mode"] == "pe":
        return {'R1':str(indir) + "/{wildcards.sample}/{wildcards.sample}_R1.qc.humanDecontaminated.fastq.gz".format(wildcards=wildcards), 'R2': str(indir) + "/{wildcards.sample}/{wildcards.sample}_R2.qc.humanDecontaminated.fastq.gz".format(wildcards=wildcards)}
    elif config["mode"] == "se":
        return {'R1':str(indir) + "/{wildcards.sample}/{wildcards.sample}_R1.qc.humanDecontaminated.fastq.gz".format(wildcards=wildcards)}
    elif config["mode"] == "rarify":
        return {'R1':str(rare_indir) + "/{wildcards.sample}/{wildcards.sample}_R1.qc.humanDecontaminated.fastq.gz".format(wildcards=wildcards)}
    else:
        raise Exception('Invalid config["mode"]')

rule diamond:
    input:
        unpack(input_reads)
    output:
        "%s/{sample}/{sample}_Diamond_{DB_name,[A-Za-z0-9]+}.tsv" %(outdir)
    params:
        N_MISMATCH = config["n_mismatch"],
        n_cores=config["n_cores"]
    resources:
        mem="8G",
        cpus=config["n_cores"]
    run:
        shell("diamond blastx --db %s -q {input.R1} -f 6 -o {output} --threads {params.n_cores}" %(config["DB_path"]))


rule diamond_all:
    input:
        expand("%s/{sample}/{sample}_Diamond_{DB_name}.tsv" %(outdir), sample = SAMPLES, DB_name=config["DB_name"])


rule diamond_combine:
    input:
        expand("%s/{sample}/{sample}_Diamond_{DB_name}.tsv" %(outdir), sample = SAMPLES, DB_name=config["DB_name"])
    output:
        "%s/%s_Diamond_%s_%s_%s_%s.csv" %(outdir, config["expt_name"], config["DB_name"], config["diamond_method"],config["min_identity_diamond"], config["tax_info"])
    params:
        diamond_method=config["diamond_method"],
        min_identity_diamond=config["min_identity_diamond"],
        tax_info=config["tax_info"],
        infile_dict=config["species_map"],
    resources:
        time="24:00:00",
        mem="8G"
    shell:
        "python code/combine_diamond.py {input} -o {output} -m {params.diamond_method} --min_identity {params.min_identity_diamond} --tax_info {params.tax_info} --infile_dict {params.infile_dict}"

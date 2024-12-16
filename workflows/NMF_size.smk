#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     NMF size Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing  #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

import pandas as pd
import re
import numpy as np
import glob
import os
import itertools as it

# Values directly from configuration file
samples = config["samples"]
COHORTS = config["cohorts"]
threads = config["threads"]
FRAG_DISTROS = config["frag_distro"]
frag_threads = config["threads"]
genome_fasta = config["genome_fasta"]
frag_repo = config["frag_repo"]
PE1 = config["PE1"]
PE2 = config["PE2"]
fasta_base = config["fasta_base"]
lib_path = config["raw_fastq"]
frag_container = config["frag_repo"]
normal_control = config.get("normal_control", "healthy")  # Default to "healthy" if not defined

# Directory values derived from data_dir in configuration YAML
data_dir             = config["data_dir"]
frag                 = data_dir + "/analysis/cfdna_wgs/frag"
frag_bams            = config["frag_bams"]
frag_fastqs          = config["raw_fastq"]
frag_frag            = data_dir + "/analysis/cfdna_wgs/frag/delfi"
frag_frag_beds       = frag_frag  + "/beds"
frag_frag_counts     = frag_frag + "/counts"
frag_frag_gc_distros = frag_frag  + "/distros"
qcdir                = data_dir + "/analysis/cfdna_wgs/qc"
benchdir             = data_dir + "/benchmark"
logdir               = data_dir + "/logs"
refdir               = config["ref_dir"]
bwa_dir              = config["ref_dir"]
frag_scriptdir       = config["frag_repo"] +  "/scripts"
samplelist_dir       = config["data_dir"] + "/samplelist/features"

# NMF specific directories/inputs
nmf_raw_hist_dir     = config["data_dir"] + "/analysis/cfdna_wgs/frag/nmf/raw"
nmf_agg_hist_dir     = config["data_dir"] + "/analysis/cfdna_wgs/frag/nmf/agg"
nmf_dir     = config["data_dir"] + "/analysis/cfdna_wgs/frag/nmf"
nmf_probs = [f"{nmf_dir}/probs/{c1}_v_{c2}.n={{nc}}.csv" for (c1,c2) in it.combinations(COHORTS,2)]
nmf_scores = [f"{nmf_dir}/scores/{c1}_v_{c2}.n={{nc}}.csv" for (c1,c2) in it.combinations(COHORTS,2)]
nmf_comps = [f"{nmf_dir}/model/{c1}_v_{c2}.n={{nc}}.pkl" for (c1,c2) in it.combinations(COHORTS,2)]
max_components = config["max_nmf_components"]
nc = range(2,config["max_nmf_components"]+1)

# Read samples
frag_libs = pd.read_table(config["samples"])

# Ensure readable fastqs
readable = []
#for x in lib_path + "/" + frag_libs["r1_basename"]:
for x in frag_libs.file:
    readable.append(os.access(x, os.R_OK))
frag_libs['readable']=readable
frag_libs = frag_libs[frag_libs.readable == True]

# Make the dictionary
FRAG_LIBS = frag_libs["library"].tolist()
frag_libs_file_indict = frag_libs["file"].tolist()
frag_lib_dict = dict(zip(FRAG_LIBS, frag_libs_file_indict))

# Make  a list of healthy libraries
FRAG_HEALTHY_LIBRARIES = frag_libs[frag_libs['cohort'] == normal_control]['library'].tolist()

#Define cohort comparisons (one-versus-one)
COHORT_1 = []
COHORT_2 = []
for c1, c2 in it.combinations(COHORTS,2):
    COHORT_1.append(c1)
    COHORT_2.append(c2)
    print(c1,c2)

rule all:
  input: 
    # Raw histograms
    expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=frag_libs["library"]),
    expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=frag_libs["library"]),
    # Aggregated histograms
    nmf_agg_hist_dir + "/count_histogram.csv",
    nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
    # NMF model outputs
    expand(nmf_dir + "/probs/{c1}_v_{c2}.n={nc}.csv", c1=COHORT_1, c2=COHORT_2, nc=range(2, config["max_nmf_components"] + 1)),
    expand(nmf_dir + "/model/{c1}_v_{c2}.n={nc}.pkl", c1=COHORT_1, c2=COHORT_2, nc=range(2, config["max_nmf_components"] + 1)),
    expand(nmf_dir + "/coeffs/{c1}_v_{c2}.n={nc}.csv", c1=COHORT_1, c2=COHORT_2, nc=range(2, config["max_nmf_components"] + 1)),
    # aucs
    nmf_dir + "/aucs/aucs.csv",
    nmf_dir + "/aucs/best_med_aucs.csv",
    expand(nmf_scores, nc = range(2,max_components+1)),
rule raw_histograms:
    resources:
        time="3:00:00"
    input:
        bam = frag_bams + "/{lib}_filt.bam",
    output:
        metrics   = nmf_raw_hist_dir + "/{lib}_metrics.txt",
        histogram = nmf_raw_hist_dir + "/{lib}_histogram.pdf",
    log:
        logdir + "/nmf/raw_histograms/{lib}.log",
    params:
        out_dir = nmf_raw_hist_dir,
        log_dir = logdir + "/nmf/raw_histograms",
        width   = config["histogram_width"],
        picardpath = config["picardpath"],
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        module load picard
        module load R

        java -jar {params.picardpath} CollectInsertSizeMetrics \
        -I {input.bam} \
        -O {output.metrics} \
        -H {output.histogram} \
        -W {params.width} \
        &> {log}
        """
rule agg_histograms:
    resources:
        time="0:90:00"
    input:
        expand(nmf_raw_hist_dir + "/{lib}_{out}",
                lib=frag_libs["library"],
                out=["metrics.txt","histogram.pdf"],
                ),
    output:
        counts = nmf_agg_hist_dir + "/count_histogram.csv",
        freqs  = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
    log:
        logdir + "/nmf/agg_histograms.log",
    params:
        hist_dir = nmf_raw_hist_dir,
        out_dir  = nmf_agg_hist_dir,
        logs_dir  = logdir + "/nmf",
        script   = frag_scriptdir + "/agg_histograms.py",
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.logs_dir}

        python3 {params.script} \
        --hist_dir {params.hist_dir} \
        --output {params.out_dir} \
        --verbose \
        &> {log}
        """
rule nmf:
    resources:
        time="0:20:00",
        mem_gb="16g",
        cpus="8",
    input:
        nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
    output:
        probs = nmf_dir + "/probs/{c1}_v_{c2}.n={nc}.csv",
        coeffs = nmf_dir + "/coeffs/{c1}_v_{c2}.n={nc}.csv",
    log:
        logdir + "/nmf/model/{c1}_v_{c2}.n={nc}.log"
    params:
        out_dir = nmf_dir,
        log_dir = logdir + "/nmf/model",
        script  = frag_scriptdir + "/nmf_model.py",
        samplelist = samples,
        #nc = "2-"(config["max_nmf_components"]+1),
        nc = lambda wildcards: wildcards.nc,
        c1 = lambda wildcards: wildcards.c1,
        c2 = lambda wildcards: wildcards.c2,
        model = lambda wildcards: nmf_dir + f"/model/{wildcards.c1}_v_{wildcards.c2}.n={wildcards.nc}.pkl",
        name = lambda wildcards: f"{wildcards.c1}_v_{wildcards.c2}.n={wildcards.nc}",
    shell:
        """
        mkdir -p {params.out_dir} {params.log_dir}

        # Check for existing model
        if [ -f {params.model} ]; then
            echo "Using existing model: {params.model}"
            python3 {params.script} \
                --freq_csv {input} \
                --samplelist {params.samplelist} \
                --n_components {wildcards.nc} \
                --cohort1 {params.c1} \
                --cohort2 {params.c2} \
                --out_dir {params.out_dir} \
                --existing_model {params.model} \
                --name {params.name} \
                &> {log}
        else
            echo "Model does not exist, creating new model"
            python3 {params.script} \
                --freq_csv {input} \
                --samplelist {params.samplelist} \
                --n_components {params.nc} \
                --cohort1 {params.c1} \
                --cohort2 {params.c2} \
                --out_dir {params.out_dir} \
                --name {params.name} \
                &> {log}
        fi
        """
rule nmf_aucs:
    resources:
        time="0:05:00",
    input:
        expand(nmf_probs,
                nc = range(2,max_components+1)
        ),
        samplelist = samples,
    output:
        nmf_dir + "/aucs/aucs.csv",
        nmf_dir + "/aucs/best_med_aucs.csv",
        expand(nmf_scores,
                nc = range(2,max_components+1)
        ),
    log:
        logdir + "/nmf/auc/aucs.log"
    params:
        out_dir = nmf_dir,
        log_dir = logdir + "/nmf/auc",
        script  = frag_scriptdir + "/nmf_aucs.py",
        probs_dir = nmf_dir + "/probs",
        max_components = max_components,
        cohorts = ",".join(config["cohorts"]),
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --probs_dir {params.probs_dir} \
        --max_components {params.max_components} \
        --samplelist {input.samplelist} \
        --out_dir {params.out_dir} \
        --cohorts {params.cohorts} \
        &> {log}
        """
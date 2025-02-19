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
bwa_dir             = config["ref_dir"]
COHORTS             = config["cohorts"]
data_dir            = config["data_dir"]
fasta_base          = config["fasta_base"]
FRAG_DISTROS        = config["frag_distro"]
frag_repo           = config["frag_repo"]
frag_threads        = config["threads"]
genome_fasta        = config["genome_fasta"]
metadata_file       = config["metadata_file"]
normal_control      = config.get("normal_control", "healthy")  # Default to "healthy" if not defined
PE1                 = config.get("PE1", "need to define PE1")
PE2                 = config.get("PE2", "need to define PE2")
refdir              = config["ref_dir"]
samples             = config["samples"]
threads             = config["threads"]

# Directory values derived from configuration YAML
benchdir            = data_dir + "/benchmark"
frag                = data_dir + "/analysis/cfdna_wgs"
frag_bams           = data_dir + "/bams"
frag_fastqs         = data_dir + "/fastq"
frag_scriptdir      = frag_repo +  "/scripts"
logdir              = data_dir + "/logs/NMF_size"
qcdir               = data_dir + "/qc"

# NMF specific inputs
max_components      = config["max_nmf_components"]
nc                  = range(2,config["max_nmf_components"]+1)

# NMF specific directories
nmf_dir             =  frag + "/nmf"

nmf_agg_hist_dir    = nmf_dir + "/agg"
nmf_comps           = [f"{nmf_dir}/model/{c1}_v_{c2}.n={{nc}}.pkl" for (c1,c2) in it.combinations(COHORTS,2)]
nmf_probs           = [f"{nmf_dir}/probs/{c1}_v_{c2}.n={{nc}}.csv" for (c1,c2) in it.combinations(COHORTS,2)]
nmf_raw_hist_dir    = nmf_dir + "/raw"
nmf_scores          = [f"{nmf_dir}/scores/{c1}_v_{c2}.n={{nc}}.csv" for (c1,c2) in it.combinations(COHORTS,2)]

####------ Read samples ------####
frag_libs = pd.read_table(config["samples"])

# Ensure readable fastqs
readable = []
for x in frag_libs.file:
    readable.append(os.access(x, os.R_OK))
frag_libs['readable']=readable
frag_libs = frag_libs[frag_libs.readable == True]

# Make the dictionary
FRAG_LIBS = frag_libs["library"].tolist()
print("FRAG_LIBS:", FRAG_LIBS)
frag_libs_file_indict = frag_libs["file"].tolist()
frag_lib_dict = dict(zip(FRAG_LIBS, frag_libs_file_indict))

#Check that all BAM files exist
print("Checking BAM files:")
for lib in FRAG_LIBS:
    bam_path = f"{frag_bams}/{lib}_filt.bam"
    print(f"  {bam_path} exists: {os.path.exists(bam_path)}")

# Make  a list of healthy libraries
FRAG_HEALTHY_LIBRARIES = frag_libs[frag_libs['cohort'] == normal_control]['library'].tolist()

#Define cohort comparisons (one-versus-one)
COHORT_1 = []
COHORT_2 = []
for c1, c2 in it.combinations(COHORTS,2):
    COHORT_1.append(c1)
    COHORT_2.append(c2)
    print(c1,c2)

def get_cohort_pairs():
    pairs = []
    for c1, c2 in it.combinations(COHORTS, 2):
        pairs.append((c1, c2))
    return pairs

rule all:
    input:
        # Raw histograms
        expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
        expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS),
        # Aggregated histograms
        nmf_agg_hist_dir + "/count_histogram.csv",
        nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
        # NMF outputs for each cohort pair and component number
        expand(nmf_dir + "/probs/{c1}_v_{c2}.n={nc}.csv",
            c1=COHORT_1, c2=COHORT_2, nc=range(2, config["max_nmf_components"] + 1)),
        expand(nmf_dir + "/model/{c1}_v_{c2}.n={nc}.pkl",
            c1=COHORT_1, c2=COHORT_2, nc=range(2, config["max_nmf_components"] + 1)),
        expand(nmf_dir + "/coeffs/{c1}_v_{c2}.n={nc}.csv",
            c1=COHORT_1, c2=COHORT_2, nc=range(2, config["max_nmf_components"] + 1))

rule raw_histograms:
    resources:
        time="5:00:00"
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
        time="0:120:00"
    input:
        metrics = expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
        plots = expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS)
    output:
        counts = nmf_agg_hist_dir + "/count_histogram.csv",
        freqs = nmf_agg_hist_dir + "/freq_histogram.30_700.csv"
    log:
        logdir + "/nmf/agg_histograms.log"
    params:
        hist_dir = nmf_raw_hist_dir,
        out_dir = nmf_agg_hist_dir,
        logs_dir = logdir + "/nmf",
        script = frag_scriptdir + "/agg_histograms.py"
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
        time="0:30:00",
        mem_gb="16g",
        cpus="8",
    input:
        freq_csv = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
        counts = nmf_agg_hist_dir + "/count_histogram.csv",

        raw_metrics = expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
        raw_plots = expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS)
    output:
        probs = nmf_dir + "/probs/{c1}_v_{c2}.n={nc}.csv",
        coeffs = nmf_dir + "/coeffs/{c1}_v_{c2}.n={nc}.csv",
        model = nmf_dir + "/model/{c1}_v_{c2}.n={nc}.pkl"
    log:
        logdir + "/nmf/model/{c1}_v_{c2}.n={nc}.log"
    params:
        out_dir = nmf_dir,
        log_dir = logdir + "/nmf/model",
        script = frag_scriptdir + "/nmf_model.py",
        samplelist = samples,
        existing_model = config["existing_nmf_model"]
    shell:
        """
        mkdir -p {params.out_dir}/{{model,probs,coeffs}} {params.log_dir}

        python3 {params.script} \
            --freq_csv {input.freq_csv} \
            --samplelist {params.samplelist} \
            --n_components {wildcards.nc} \
            --cohort1 {wildcards.c1} \
            --cohort2 {wildcards.c2} \
            --out_dir {params.out_dir} \
            --existing_model {params.existing_model} \
            --name {wildcards.c1}_v_{wildcards.c2}.n={wildcards.nc} \
            &> {log}
        """
rule nmf_aucs:
    resources:
        time="0:25:00",
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
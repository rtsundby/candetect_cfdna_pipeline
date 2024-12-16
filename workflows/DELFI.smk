#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
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
cfdna_wgs_container =  config["frag_repo"]
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
frag_scriptdir = config["frag_repo"] +  "/scripts"
samplelist_dir    = config["data_dir"] + "/samplelist/features"

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
#filt_bam_to_frag_bed
    expand(frag_frag_beds + "/{library}_filt.bed",
        library = frag_libs["library"],
    ),
#gc_distro
    expand(frag_frag_gc_distros + "/{library}_gc_distro.csv",
        library = frag_libs["library"],
    ),
#healthy_gc
    frag_frag_gc_distros + "/healthy_med.rds",
#frag_gc_sample
    expand(frag_frag_beds + "/{library}_sampled_frag.bed",
        library = frag_libs["library"],
    ),
#frag_sum
    expand(frag_frag_beds + "/{library}_norm_short.bed",
        library = frag_libs["library"],
    ),
    expand(frag_frag_beds + "/{library}_norm_long.bed",
        library = frag_libs["library"],
    ),
#frag_window_count
    expand(frag_frag_counts + "/{library}_cnt_short.tmp",
        library = frag_libs["library"],
    ),
    expand(frag_frag_counts + "/{library}_cnt_long.tmp",
        library = frag_libs["library"],
    ),
#frag_count_merge
    frag_frag + "/frag_counts.tsv",
#unit_cent_sd
    frag_frag + "/ratios.tsv",
# delfi armz (raw)
    frag_frag + "/delfi_armz.raw.csv",
# delfi_features
    frag_frag + "/delfi_armz.clean.csv",
    frag_frag + "/delfi_ratio.norm_rows.clean.csv",
# delfi_model
    expand(frag_frag + "/probs/{c1}_v_{c2}.csv", 
    zip,
        c1=COHORT_1,
        c2=COHORT_2
    ),
    expand(frag_frag + "/coeffs/{c1}_v_{c2}.csv", 
    zip,
        c1=COHORT_1,
        c2=COHORT_2
    ),
# aucs
    frag_frag + "/aucs/aucs.csv",
    expand(frag_frag + "/scores/{c1}_v_{c2}.csv",
    zip,
        c1 = COHORT_1,
        c2 = COHORT_2
    ),
#prob plots
    logdir + "/delfi/plot_probs.log",
# coefficient heatmaps
    logdir + "/delfi/delfi_plot_coeff_heatmaps.log",

rule filt_bam_to_frag_bed:
    benchmark: benchdir + "/{library}_filt_bam_to_frag_bed.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/delfi/bam_to_bed/{library}_filt_bam_to_frag_bed.log",
    output: frag_frag_beds + "/{library}_filt.bed",
    params:
        fragdir = frag,
        frag_fragdir = frag_frag,
        out_dir = frag_frag_beds,
        fasta = genome_fasta,
        script = frag_scriptdir + "/filt_bam_to_frag_bed.sh",
        threads = frag_threads,
        log_dir = logdir + "/delfi/bam_to_bed",
    shell:
        """
        mkdir -p {params.log_dir} 
        mkdir -p {params.fragdir} \
        mkdir -p {params.frag_fragdir} \
        mkdir -p {params.out_dir} \

        {params.script} \
        {input} \
        {params.fasta} \
        {params.threads} \
        {output}
        """

# Make GC distributions
rule gc_distro:
    benchmark: benchdir + "/{library}_frag_gc_distro.benchmark.txt",
    input: frag_frag_beds + "/{library}_filt.bed",
    log: logdir + "/delfi/gc_distro/{library}_frag_gc_distro.log",
    output: frag_frag_gc_distros + "/{library}_gc_distro.csv",
    params:
        script = frag_scriptdir + "/gc_distro.R",
        out_dir = frag_frag_gc_distros,
        log_dir = logdir + "/delfi/gc_distro",
    shell:
        """
        mkdir -p {params.log_dir} 
        mkdir -p {params.out_dir} \
        module load R
        Rscript {params.script} \
        {input} \
        {output} \
        > {log} 2>&1
        """

# Make healthy GC distributions summary file
rule healthy_gc:
    benchmark: benchdir + "/frag_healthy_gc.benchmark.txt",
    input: expand(frag_frag_gc_distros + "/{library}_gc_distro.csv", library = FRAG_HEALTHY_LIBRARIES),
    log: logdir + "/delfi/healthy_gc/frag_healthy_gc.log",
    output: frag_frag_gc_distros + "/healthy_med.rds",
    params:
        distro_dir = frag_frag_gc_distros,
        log_dir = logdir + "/delfi/healthy_gc",
        script = frag_scriptdir + "/healthy_gc.R",
    shell:
        """
        mkdir -p {params.log_dir} 
        module load R
        Rscript {params.script} \
        {params.distro_dir} \
        "{input}" \
        {output} > {log} 2>&1
        """

# Sample fragments by healthy GC proportions
rule frag_gc_sample:
    benchmark: benchdir + "/{library}_frag_gc_sample.benchmark.txt",
    input:
        frag_bed = frag_frag_beds + "/{library}_filt.bed",
        healthy_med = frag_frag_gc_distros + "/healthy_med.rds",
    log: logdir + "/delfi/gc_sample/{library}_frag_gc_sample.log",
    output: frag_frag_beds + "/{library}_sampled_frag.bed",
    params:
        log_dir = logdir + "/delfi/gc_sample",
        out_dir = frag_frag_beds,
        script = frag_scriptdir + "/gc_sample.R",
    shell:
        """
        mkdir -p {params.log_dir} 
        mkdir -p {params.out_dir} \

        module load R
        Rscript {params.script} \
        {input.healthy_med} \
        {input.frag_bed} \
        {output} > {log} 2>&1
        """
# Sum fragments in short and long length groups
rule frag_sum:
    benchmark: benchdir + "/{library}_frag_sum.benchmark.txt",
    input: frag_frag_beds + "/{library}_sampled_frag.bed",
    log: logdir + "/delfi/frag_sum/{library}_frag_frag_window_sum.log",
    output:
        short = frag_frag_beds + "/{library}_norm_short.bed",
        long =  frag_frag_beds + "/{library}_norm_long.bed",
    params:
        log_dir = logdir + "/delfi/frag_sum",
        script = frag_scriptdir + "/frag_window_sum.sh",
        threads = frag_threads,
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {output.short} {output.long} &> {log}
        """

# Count short and long fragments intersecting kept genomic windows
rule frag_window_count:
    benchmark: benchdir + "/{library}_frag_frag_window_int.benchmark.txt",
    input:
        short = frag_frag_beds + "/{library}_norm_short.bed",
        long = frag_frag_beds + "/{library}_norm_long.bed",
        matbed = refdir + "/keep_5mb.bed",
    log: logdir + "/delfi/frag_window_count/{library}_frag_frag_window_int.log",
    output:
        short = frag_frag_counts + "/{library}_cnt_short.tmp",
        long = frag_frag_counts + "/{library}_cnt_long.tmp",
    params:
        script = frag_scriptdir + "/frag_window_int.sh",
        threads = threads,
        out_dir = frag_frag_counts,
        log_dir = logdir + "/delfi/frag_window_count",
    shell:
        """
        mkdir -p {params.out_dir} 
        mkdir -p {params.log_dir}
        {params.script} \
        {input.short} \
        {input.matbed} \
        {output.short}
        {params.script} \
        {input.long} \
        {input.matbed} \
        {output.long}
        """

# Merge short and long fragment counts by genomic poistion for all libraries
rule frag_count_merge:
    benchmark: benchdir + "/frag_count_merge.benchmark.txt",
    input: expand(frag_frag_counts + "/{library}_cnt_{length}.tmp",  library = FRAG_LIBS, length = ["short","long"]),
    log: logdir + "/delfi/frag_count_merge/frag_count_merge.log",
    output:  frag_frag + "/frag_counts.tsv",
    params:
        counts_dir = frag_frag + "/counts",
        script = frag_scriptdir + "/count_merge.sh",
        threads = frag_threads,
        log_dir = logdir + "/delfi/frag_count_merge",
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {params.counts_dir} \
        {output} &> {log}
        """
rule unit_cent_sd:
    benchmark: benchdir + "/unit_cent_sd.benchmark.txt",
    input: frag_frag + "/frag_counts.tsv",
    log: logdir + "/delfi/unit_cent_sd/unit_cent_sd.log",
    output: frag_frag + "/ratios.tsv",
    params:
        script = frag_scriptdir + "/make_ratios.R",
        log_dir = logdir + "/delfi/unit_cent_sd",
    shell:
        """
        mkdir -p {params.log_dir}
        module load R
        Rscript {params.script} \
        {input} {output} > {log} 2>&1
        """
rule delfi_armz_raw:
    resources:
        time="0:60:00",
        mem_gb="4g",
    input:
        frag_frag + "/frag_counts.tsv",
    output:
        frag_frag + "/delfi_armz.raw.csv",
    log:
        logdir + "/delfi/armz_raw/delfi_armz_raw.log",
    params:
        log_dir = logdir + "/delfi/armz_raw",
        script = frag_scriptdir + "/delfi_armz.R",
        cytoband = config["cytoband"],
        sample_list = samples,
        control_term = normal_control,
    shell:
        """
        mkdir -p {params.log_dir}
        module load R
        
        Rscript {params.script} \
        {params.cytoband} \
        {params.sample_list} \
        {params.control_term} \
        {input} \
        {output} \
        &> {log}
        """
rule delfi_features:
    resources:
        time="0:60:00",
    input:
        arm_z = frag_frag +  "/delfi_armz.raw.csv",
        ratio = frag_frag + "/ratios.tsv",
    output:
        arm_z = frag_frag + "/delfi_armz.clean.csv",
        ratio = frag_frag + "/delfi_ratio.norm_rows.clean.csv",
    log:
        logdir + "/delfi/features/delfi_features.log",
    params:
        script  = frag_scriptdir + "/delfi_features.py",
        out_dir = frag_frag,
        log_dir = logdir + "/delfi/features",
    shell:
        """
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --armz_csv {input.arm_z} \
        --ratio_tsv {input.ratio} \
        --out_dir {params.out_dir} \
        &> {log}
        """      
rule delfi_model:
    resources:
        time="0:10:00",
        mem_gb="2g",
        cpus="8",
    input:
       frag_frag + "/delfi_ratio.norm_rows.clean.csv",
    output:
        probs = frag_frag + "/probs/{c1}_v_{c2}.csv",
        coeffs = frag_frag + "/coeffs/{c1}_v_{c2}.csv",
    log:
        logdir + "/delfi/model/delfi_model_{c1}_v_{c2}.log"
    params:
        arm_z = frag_frag + "/delfi_armz.clean.csv",
        sample_list = samples,
        model = lambda wildcards: frag_frag + f"/model/{wildcards.c1}_v_{wildcards.c2}.pkl",
        out_dir = frag_frag,
        out_dirprob = frag_frag + "/probs",
        out_dircoeffs = frag_frag + "/coeffs",
        script = frag_scriptdir + "/delfi_model.py",
        c1 = lambda wildcards: wildcards.c1,
        c2 = lambda wildcards: wildcards.c2,
        log_dir = logdir + "/delfi",
        use_existing_model = lambda wildcards: frag_frag + f"/model/{wildcards.c1}_v_{wildcards.c2}.pkl" if os.path.exists(frag_frag + f"/model/{wildcards.c1}_v_{wildcards.c2}.pkl") else None,
        name = "{c1}_v_{c2}",
        log_dir = logdir + "/delfi/model",
    shell:
        """
        mkdir -p {params.out_dir} {params.out_dirprob} {params.out_dircoeffs} {params.log_dir}

        # Use conditional for existing model
        if [ -f {params.model} ]; then
            echo "Using existing model: {params.model}"
            python3 {params.script} \
                --armz_csv {params.arm_z} \
                --ratio_csv {input} \
                --samplelist {params.sample_list} \
                --cohort1 {params.c1} \
                --cohort2 {params.c2} \
                --out_dir {params.out_dir} \
                --name {params.name} \
                --existing_model {params.model} \
                &> {log}
        else
            echo "Model does not exist, creating new model"
            python3 {params.script} \
                --armz_csv {params.arm_z} \
                --ratio_csv {input} \
                --samplelist {params.sample_list} \
                --cohort1 {params.c1} \
                --cohort2 {params.c2} \
                --out_dir {params.out_dir} \
                --name {params.name} \
                &> {log}
        fi
        """
rule delfi_aucs:
    resources:
        time="0:05:00",
    input:
        expand(frag_frag + "/probs/{c1}_v_{c2}.csv",
        zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        sample_list = samples,
    output:
        frag_frag + "/aucs/aucs.csv",
        expand(frag_frag + "/scores/{c1}_v_{c2}.csv",
                zip,
                c1 = COHORT_1,
                c2 = COHORT_2
        ),
    log:
        logdir + "/delfi/auc/delfi_aucs.log"
    params:
        out_dir = frag_frag,
        out_dir_auc = frag_frag + "/aucs",
        out_dir_scores = frag_frag + "/scores",
        script = frag_scriptdir + "/aucs.py",
        probs_dir = frag_frag + "/probs",
        log_dir = logdir + "/delfi/auc",
    shell:
        """
        mkdir -p {params.out_dir_auc} {params.out_dir_scores} {params.log_dir} \

        python3 {params.script} \
        --probs_dir {params.probs_dir} \
        --samplelist {input.sample_list} \
        --out_dir {params.out_dir} \
        &> {log}
        """
rule plot_probs:
    resources:
        time="0:20:00", 
    input:
        expand(frag_frag + "/probs/{c1}_v_{c2}.csv",
        zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        sample_list = samples,
    log:
        logdir + "/delfi/plot_probs/plot_probs.log"
    params:
        out_dir = frag_frag,
        log_dir = logdir + "/delfi/plot_probs",
        script  = frag_scriptdir + "/plot_probs.py",
        probs_dir = frag_frag + "/probs",
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --probs_dir {params.probs_dir} \
        --samplelist {input.sample_list} \
        --out_dir {params.out_dir} \
        &> {log}
        """
rule plot_coeff_heatmaps:
    resources:
        time="0:20:00",
    input:
        arm_z = frag_frag + "/delfi_armz.clean.csv",
    log:
        logdir + "/delfi/coeff_heatmaps/delfi_plot_coeff_heatmaps.log",
    params:
        log_dir = logdir + "/delfi/coeff_heatmaps",
        script  = frag_scriptdir + "/delfi_plot_coefficient_heatmaps.py",
        data_dir = frag_frag,
    shell:
        """
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --data_dir {params.data_dir} \
        --armz_csv {input.arm_z} \
        &> {log}
        """
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

# Directory values derived from data_dir in configuration YAML
benchdir            = data_dir + "/benchmark"
frag                = data_dir + "/analysis/cfdna_wgs"
frag_bams           = data_dir + "/bams"
frag_fastqs         = data_dir + "/fastq"
frag_scriptdir      = frag_repo +  "/scripts"
logdir              = data_dir + "/logs/bin_wise_fragmentomics"
qcdir               = data_dir + "/qc"

# Bin-wise fragmentomic specific directories
bin_wise_frag            = frag + "/bin_wise_fragmentomics"

bin_wise_beds      = bin_wise_frag + "/beds"
bin_wise_counts     = bin_wise_frag + "/counts"
bin_wise_gc_distros = bin_wise_frag + "/distros"

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

# Mode 3: apply existing model + existing cutpoints (no new AUC)
APPLY_BINWISE_CUTPOINTS = config.get("existing_binwise_scores_dir", "no") != "no"

# Shared preprocessing and model outputs (all modes)
_binwise_base_outputs = [
    expand(bin_wise_beds + "/{library}_filt.bed", library=frag_libs["library"]),
    expand(bin_wise_gc_distros + "/{library}_gc_distro.csv", library=frag_libs["library"]),
    bin_wise_gc_distros + "/healthy_med.rds",
    expand(bin_wise_beds + "/{library}_sampled_frag.bed", library=frag_libs["library"]),
    expand(bin_wise_beds + "/{library}_norm_short.bed", library=frag_libs["library"]),
    expand(bin_wise_beds + "/{library}_norm_long.bed", library=frag_libs["library"]),
    expand(bin_wise_counts + "/{library}_cnt_short.tmp", library=frag_libs["library"]),
    expand(bin_wise_counts + "/{library}_cnt_long.tmp", library=frag_libs["library"]),
    bin_wise_frag + "/frag_counts.tsv",
    bin_wise_frag + "/ratios.tsv",
    bin_wise_frag + "/binwise_armz.raw.csv",
    bin_wise_frag + "/binwise_armz.clean.csv",
    bin_wise_frag + "/binwise_ratio.norm_rows.clean.csv",
    expand(bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv", zip, c1=COHORT_1, c2=COHORT_2),
    expand(bin_wise_frag + "/model/coeffs/{c1}_v_{c2}.csv", zip, c1=COHORT_1, c2=COHORT_2),
]
# Mode 3: scores only via cutpoint application
if APPLY_BINWISE_CUTPOINTS:
    _binwise_analysis_outputs = [
        expand(bin_wise_frag + "/scores/{c1}_v_{c2}.csv", zip, c1=COHORT_1, c2=COHORT_2),
    ]
# Modes 1/2: full AUC + plots
else:
    _binwise_analysis_outputs = [
        bin_wise_frag + "/aucs/aucs/aucs.csv",
        expand(bin_wise_frag + "/scores/{c1}_v_{c2}.csv", zip, c1=COHORT_1, c2=COHORT_2),
        expand(bin_wise_frag + "/plots/prob_per_repeat/{c1}_v_{c2}.pdf", zip, c1=COHORT_1, c2=COHORT_2),
        expand(bin_wise_frag + "/plots/prob_per_sample/{c1}_v_{c2}.pdf", zip, c1=COHORT_1, c2=COHORT_2),
        bin_wise_frag + "/binwise/plots/coeff_heatmap/binwise_armz.clean.pdf",
    ]

rule all:
  input:
    _binwise_base_outputs,
    _binwise_analysis_outputs,

rule filt_bam_to_frag_bed:
    resources:
        time="72:00:00",
        mem_mb="100000",
        cpus="30",
    benchmark: benchdir + "/{library}_filt_bam_to_frag_bed.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/bam_to_bed/{library}_filt_bam_to_frag_bed.log",
    output: bin_wise_beds + "/{library}_filt.bed",
    params:
        fragdir = frag,
        bin_wise_fragdir = bin_wise_frag,
        out_dir = bin_wise_beds,
        fasta = genome_fasta,
        script = frag_scriptdir + "/filt_bam_to_frag_bed.sh",
        threads = frag_threads,
        log_dir = logdir + "/bam_to_bed",
    shell:
        """
        mkdir -p {params.log_dir} 
        mkdir -p {params.fragdir} 
        mkdir -p {params.bin_wise_fragdir} 
        mkdir -p {params.out_dir} 

        {params.script} \
        {input} \
        {params.fasta} \
        {params.threads} \
        {output}
        """

# Make GC distributions
rule gc_distro:
    resources:
        time="12:00:00",
        mem_mb="64000",
    benchmark: benchdir + "/{library}_frag_gc_distro.benchmark.txt",
    input: bin_wise_beds + "/{library}_filt.bed",
    log: logdir + "/gc_distro/{library}_frag_gc_distro.log",
    output: bin_wise_gc_distros + "/{library}_gc_distro.csv",
    params:
        script = frag_scriptdir + "/gc_distro.R",
        out_dir = bin_wise_gc_distros,
        log_dir = logdir + "/gc_distro",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir}

        module load R
        Rscript {params.script} \
        {input} \
        {output} > {log} 2>&1
        """

# Make healthy GC distributions summary file
rule healthy_gc:
    resources:
        time="12:00:00",
        mem_mb="32000",
    benchmark: benchdir + "/frag_healthy_gc.benchmark.txt",
    input: expand(bin_wise_gc_distros + "/{library}_gc_distro.csv", library = FRAG_HEALTHY_LIBRARIES),
    log: logdir + "/healthy_gc/frag_healthy_gc.log",
    output: bin_wise_gc_distros + "/healthy_med.rds",
    params:
        distro_dir = bin_wise_gc_distros,
        log_dir = logdir + "/healthy_gc",
        script = frag_scriptdir + "/healthy_gc.R",
    shell:
        """
        mkdir -p {params.log_dir} 
        module load R && \
        Rscript {params.script} \
        {params.distro_dir} \
        "{input}" \
        {output} > {log} 2>&1
        """

# Sample fragments by healthy GC proportions
rule frag_gc_sample:
    resources:
        time="12:00:00",
        mem_mb="32000",
    benchmark: benchdir + "/{library}_frag_gc_sample.benchmark.txt",
    input:
        frag_bed = bin_wise_beds+ "/{library}_filt.bed",
        healthy_med = bin_wise_gc_distros + "/healthy_med.rds",
    log: logdir + "/gc_sample/{library}_frag_gc_sample.log",
    output: bin_wise_beds+ "/{library}_sampled_frag.bed",
    params:
        log_dir = logdir + "/gc_sample",
        out_dir = bin_wise_beds,
        script = frag_scriptdir + "/gc_sample.R",
    shell:
        """
        mkdir -p {params.log_dir} 
        mkdir -p {params.out_dir} 

        module load R
        Rscript {params.script} \
        {input.healthy_med} \
        {input.frag_bed} \
        {output} > {log} 2>&1
        """
# Sum fragments in short and long length groups
rule frag_sum:
    resources:
        time="8:00:00",
        mem_mb="64000",
        cpus="16",
    benchmark: benchdir + "/{library}_frag_sum.benchmark.txt",
    input: bin_wise_beds+ "/{library}_sampled_frag.bed",
    log: logdir + "/frag_sum/{library}_bin_wise_frag_window_sum.log",
    output:
        short = bin_wise_beds+ "/{library}_norm_short.bed",
        long =  bin_wise_beds+ "/{library}_norm_long.bed",
    params:
        log_dir = logdir + "/frag_sum",
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
    resources:
        time="12:00:00",
        mem_mb="200000",
        cpus="40",
    benchmark: benchdir + "/{library}_bin_wise_frag_window_int.benchmark.txt",
    input:
        short = bin_wise_beds+ "/{library}_norm_short.bed",
        long = bin_wise_beds+ "/{library}_norm_long.bed",
        matbed = refdir + "/keep_5mb.bed",
        genome = genome_fasta + ".fai",
    log: logdir + "/frag_window_count/{library}_bin_wise_frag_window_int.log",
    output:
        short = bin_wise_counts + "/{library}_cnt_short.tmp",
        long = bin_wise_counts + "/{library}_cnt_long.tmp",
    params:
        script = frag_scriptdir + "/frag_window_int.sh",
        threads = threads,
        out_dir = bin_wise_counts,
        log_dir = logdir + "/frag_window_count",
    shell:
        """
        # Create output directories
        mkdir -p {params.out_dir} 
        mkdir -p {params.log_dir}

        # Log start of processing
        echo "Starting processing for library {wildcards.library}" >> {log} 2>&1
        echo "Input files:" >> {log} 2>&1
        echo "Short fragments: {input.short}" >> {log} 2>&1
        echo "Long fragments: {input.long}" >> {log} 2>&1
        echo "Matrix bed: {input.matbed}" >> {log} 2>&1
        echo "Genome file: {input.genome}" >> {log} 2>&1

        # Process short fragments
        echo "Processing short fragments..." >> {log} 2>&1
        {params.script} \
            {input.short} \
            {input.matbed} \
            {output.short} \
            {input.genome} >> {log} 2>&1 || {{
            echo "ERROR: Failed to process short fragments for {wildcards.library}" >> {log} 2>&1
            echo "Check the log file for details: {log}" >> {log} 2>&1
            exit 1
        }}

        # Process long fragments
        echo "Processing long fragments..." >> {log} 2>&1
        {params.script} \
            {input.long} \
            {input.matbed} \
            {output.long} \
            {input.genome} >> {log} 2>&1 || {{
            echo "ERROR: Failed to process long fragments for {wildcards.library}" >> {log} 2>&1
            echo "Check the log file for details: {log}" >> {log} 2>&1
            exit 1
        }}

        # Verify outputs
        echo "Verifying output files..." >> {log} 2>&1
        for out_file in {output.short} {output.long}; do
            if [[ ! -s $out_file ]]; then
                echo "ERROR: Output file is empty: $out_file" >> {log} 2>&1
                exit 1
            fi
            echo "Successfully created: $out_file" >> {log} 2>&1
            echo "First few lines of output:" >> {log} 2>&1
            head -n 3 $out_file >> {log} 2>&1
        done

        echo "Successfully completed processing for {wildcards.library}" >> {log} 2>&1
        """

# Merge short and long fragment counts by genomic poistion for all libraries
rule frag_count_merge:
    resources:
        time="5:30:00",
        mem_mb="200000",
    benchmark:
        benchdir + "/frag_count_merge.benchmark.txt",
    input:
        expand(bin_wise_counts + "/{library}_cnt_{length}.tmp",
            library = FRAG_LIBS,
            length = ["short","long"],
        ),
    log:
        logdir + "/frag_count_merge/frag_count_merge.log",
    output:
        bin_wise_frag + "/frag_counts.tsv",
    params:
        counts_dir = bin_wise_counts,
        script = frag_scriptdir + "/count_merge.sh",
        log_dir = logdir + "/frag_count_merge",
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {params.counts_dir} \
        {output} \
        &> {log}
        """
rule unit_cent_sd:
    resources:
        time="1:30:00",
        mem_mb="30000",
    benchmark: benchdir + "/unit_cent_sd.benchmark.txt",
    input: bin_wise_frag + "/frag_counts.tsv",
    log: logdir + "/unit_cent_sd/unit_cent_sd.log",
    output: bin_wise_frag + "/ratios.tsv",
    params:
        script = frag_scriptdir + "/make_ratios.R",
        log_dir = logdir + "/unit_cent_sd",
    shell:
        """
        mkdir -p {params.log_dir}
        module load R
        Rscript {params.script} \
        {input} {output} > {log} 2>&1
        """
rule binwise_armz_raw:
    resources:
        time="0:60:00",
        mem_mb="30000",
        cpus="8",
    input:
        bin_wise_frag + "/frag_counts.tsv",
    output:
        bin_wise_frag + "/binwise_armz.raw.csv",
    log:
        logdir + "/armz_raw/binwise_armz_raw.log",
    params:
        log_dir = logdir + "/armz_raw",
        script = frag_scriptdir + "/binwise_armz.R",
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
rule binwise_features:
    resources:
        time="0:60:00",
        mem_mb="30000",
        cpus="8",
    input:
        arm_z = bin_wise_frag +  "/binwise_armz.raw.csv",
        ratio = bin_wise_frag + "/ratios.tsv",
    output:
        arm_z = bin_wise_frag + "/binwise_armz.clean.csv",
        ratio = bin_wise_frag + "/binwise_ratio.norm_rows.clean.csv",
    log:
        logdir + "/features/binwise_features.log",
    params:
        script  = frag_scriptdir + "/binwise_features.py",
        out_dir = bin_wise_frag,
        log_dir = logdir + "/features",
    shell:
        """
        mkdir -p {params.log_dir}

        python3 {params.script} \
        --armz_csv {input.arm_z} \
        --ratio_tsv {input.ratio} \
        --out_dir {params.out_dir} \
        &> {log}
        """      
rule binwise_model:
    resources:
        time="0:60:00",
        mem_mb="30000",
        cpus="8",
    input:
        ratio = bin_wise_frag + "/binwise_ratio.norm_rows.clean.csv",
    output:
        probs = bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
        coeffs = bin_wise_frag + "/model/coeffs/{c1}_v_{c2}.csv",
    log:
        logdir + "/model/binwise_model_{c1}_v_{c2}.log"
    params:
        arm_z = bin_wise_frag + "/binwise_armz.clean.csv",
        sample_list = samples,
        model = lambda wildcards: config["existing_delfi_model"] if config["existing_delfi_model"] != "no" else None,
        out_dir = bin_wise_frag + "/model",
        out_dirprob = bin_wise_frag + "/model/prob",
        out_dircoeffs = bin_wise_frag + "/model/coeffs",
        script = frag_scriptdir + "/binwise_model.py",
        c1 = lambda wildcards: wildcards.c1,
        c2 = lambda wildcards: wildcards.c2,
        name = "{c1}_v_{c2}",
        log_dir = logdir + "/model",
    shell:
        """
        mkdir -p {params.out_dir} {params.out_dirprob} {params.out_dircoeffs} {params.log_dir}

        if [ "{params.model}" != "None" ]; then
            echo "Using existing model from config: {params.model}"
            python3 {params.script} \
            --armz_csv {params.arm_z} \
            --ratio_csv {input.ratio} \
            --samplelist {params.sample_list} \
            --cohort1 {params.c1} \
            --cohort2 {params.c2} \
            --name {params.name} \
            --out_dir {params.out_dir} \
            --existing_model {params.model} \
            &> {log}
        else
            echo "Training new model"
            python3 {params.script} \
            --armz_csv {params.arm_z} \
            --ratio_csv {input.ratio} \
            --samplelist {params.sample_list} \
            --cohort1 {params.c1} \
            --cohort2 {params.c2} \
            --name {params.name} \
            --out_dir {params.out_dir} \
            &> {log}
        fi
        """
# Modes 1/2: compute AUC and plots
if not APPLY_BINWISE_CUTPOINTS:
    rule binwise_aucs:
        resources:
            time="2:00:00",
            mem_mb="30000",
            cpus="8",
        input:
            expand(bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
            zip,
                c1 = COHORT_1,
                c2 = COHORT_2
            ),
            sample_list = samples,
        output:
            bin_wise_frag + "/aucs/aucs/aucs.csv",
            expand(bin_wise_frag + "/scores/{c1}_v_{c2}.csv",
                    zip,
                    c1 = COHORT_1,
                    c2 = COHORT_2
            ),
        log:
            logdir + "/auc/binwise_aucs.log"
        params:
            out_dir = bin_wise_frag,
            out_dir_auc = bin_wise_frag + "/aucs",
            out_dir_scores = bin_wise_frag + "/scores",
            script = frag_scriptdir + "/binwise_auc.py",
            probs_dir = bin_wise_frag + "/model/probs",
            log_dir = logdir + "/auc",
            cohorts = ",".join(config["cohorts"]),
        shell:
            """
            mkdir -p {params.out_dir_auc}
            mkdir -p {params.out_dir_scores}
            mkdir -p {params.log_dir}

            python3 {params.script} \
            --probs_dir {params.probs_dir} \
            --samplelist {input.sample_list} \
            --out_dir {params.out_dir} \
            --cohorts {params.cohorts} \
            &> {log}
            """
    rule plot_probs:
        resources:
            time="0:20:00",
            mem_mb="30000",
            cpus="8",
        input:
            expand(bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
            zip,
                c1 = COHORT_1,
                c2 = COHORT_2
            ),
            sample_list = samples,
        output:
            expand(bin_wise_frag + "/plots/prob_per_repeat/{c1}_v_{c2}.pdf",
            zip,
                c1 = COHORT_1,
                c2 = COHORT_2
            ),
            expand(bin_wise_frag + "/plots/prob_per_sample/{c1}_v_{c2}.pdf",
            zip,
                c1 = COHORT_1,
                c2 = COHORT_2
            ),
        log:
            logdir + "/plot_probs/plot_probs.log"
        params:
            out_dir = bin_wise_frag,
            log_dir = logdir + "/plot_probs",
            script  = frag_scriptdir + "/plot_probs.py",
            probs_dir = bin_wise_frag + "/model/probs",
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
            mem_mb="30000",
            cpus="8",
        input:
            arm_z = bin_wise_frag + "/binwise_armz.clean.csv",
        output:
            bin_wise_frag + "/binwise/plots/coeff_heatmap/binwise_armz.clean.pdf",
        log:
            logdir + "/coeff_heatmaps/binwise_plot_coeff_heatmaps.log",
        params:
            log_dir = logdir + "/coeff_heatmaps",
            script  = frag_scriptdir + "/binwise_plot_coefficient_heatmaps.py",
            data_dir = bin_wise_frag,
        shell:
            """
            mkdir -p {params.log_dir}

            python3 {params.script} \
            --data_dir {params.data_dir} \
            --armz_csv {input.arm_z} \
            &> {log}
            """

# Mode 3: apply existing cutpoints, no new AUC
else:
    rule binwise_apply_cutpoints:
        resources:
            time="0:30:00",
            mem_mb="16000",
            cpus="4",
        input:
            probs = bin_wise_frag + "/model/probs/{c1}_v_{c2}.csv",
        output:
            scores = bin_wise_frag + "/scores/{c1}_v_{c2}.csv",
        log:
            logdir + "/apply_cutpoints/binwise_apply_cutpoints_{c1}_v_{c2}.log"
        params:
            existing_scores_dir = config.get("existing_binwise_scores_dir", "no"),
            script = frag_scriptdir + "/apply_cutpoints.py",
            samplelist = samples,
            log_dir = logdir + "/apply_cutpoints",
        shell:
            """
            mkdir -p {params.log_dir}
            mkdir -p $(dirname {output.scores})

            python3 {params.script} \
            --probs_csv {input.probs} \
            --scores_csv {params.existing_scores_dir}/scores/{wildcards.c1}_v_{wildcards.c2}.csv \
            --samplelist {params.samplelist} \
            --cohort1 {wildcards.c1} \
            --cohort2 {wildcards.c2} \
            --out_csv {output.scores} \
            &> {log}
            """
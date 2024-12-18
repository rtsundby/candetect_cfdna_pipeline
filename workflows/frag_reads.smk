##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np
import os

# Values directly from configuration file
samples = config["samples"]
threads = config["threads"]
FRAG_DISTROS = config["frag_distro"]
frag_threads = config["threads"]
genome_fasta = config["genome_fasta"]
frag_repo = config["frag_repo"]
PE1 = config.get("PE1", "need to define PE1")
print(f"PE1: {PE1}")
PE2 = config.get("PE2", "need to define PE2")
print(f"PE2: {PE2}")
fasta_base = config["fasta_base"]
lib_path = config["raw_fastq"]

# Directory values derived from data_dir in configuration YAML
data_dir             = config["data_dir"]
frag                 = data_dir + "/analysis/frag"
frag_bams            = data_dir + "/analysis/cfdna_wgs/bams"
frag_fastqs          = data_dir + "/analysis/cfdna_wgs/fastqs"
frag_frag            = data_dir + "/analysis/frag/frag"
frag_frag_beds       = data_dir + "/analysis/frag/frag/beds"
frag_frag_counts     = data_dir + "/analysis/frag/frag/counts"
frag_frag_gc_distros = data_dir + "/analysis/frag/frag/distros"
qcdir                = data_dir + "/analysis/qc"
bench_dir             = data_dir + "/benchmark"
logdir               = data_dir + "/logs/preprocessing"
refdir               = config["ref_dir"]
bwa_dir              = config["ref_dir"]
frag_scriptdir = config["frag_repo"] +  "/scripts"

#Make sure logdir exists
os.makedirs(logdir, exist_ok=True)

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
FRAG_HEALTHY_LIBRARIES = frag_libs[frag_libs['cohort'] == 'healthy']['library'].tolist()

rule all:
  input: 
#symlinks
    expand(data_dir + "/analysis/cfdna_wgs/fastqs/{library}_raw_"+ PE1 + ".fastq.gz",
            library = frag_libs["library"],
    ),
    expand(data_dir + "/analysis/cfdna_wgs/fastqs/{library}_raw_"+ PE2 + ".fastq.gz",
            library = frag_libs["library"],
    ),
#frag_index
    refdir + "/"+ fasta_base +".sa",
#frag_fastp
    expand(frag_fastqs + "/{library}_processed_"+ PE1 + ".fastq.gz",
        library = frag_libs["library"],
    ),
    expand(frag_fastqs + "/{library}_processed_"+ PE2 + ".fastq.gz",
       library = frag_libs["library"],
    ),
#frag_align
    expand(frag_bams + "/{library}_raw.bam",
        library = frag_libs["library"],
    ),
    expand(frag_bams + "/{library}_raw.bam.bai",
        library = frag_libs["library"],
    ),
#frag_dedup
    expand(frag_bams + "/{library}_dedup.bam",
        library = frag_libs["library"],
    ),
#frag_filter
    expand(frag_bams + "/{library}_filt.bam",
        library = frag_libs["library"],
    ),
#frag_fastqc
    expand(qcdir + "/{library}_{processing}_{read}_fastqc.html", 
        library = FRAG_LIBS, 
        processing = ["raw", "processed", "unpaired"], 
        read = [PE1,PE2],
        ),
    expand(qcdir + "/{library}_{processing}_{read}_fastqc.zip", 
        library = FRAG_LIBS, 
        processing = ["raw", "processed", "unpaired"], 
        read = [PE1,PE2],
        ),
#alignment_qc
    expand(qcdir + "/{library}_{processing}_samstats.txt", 
        library = FRAG_LIBS, 
        processing = ["raw","filt"],
        ),
    expand(qcdir + "/{library}_{processing}_flagstat.txt", 
        library = FRAG_LIBS, 
        processing = ["raw","filt"],
        ),
#sequencing_depth
    expand(qcdir + "/{library}_picard_depth.txt", 
        library = FRAG_LIBS,
        ),
#bampefragsize
    qcdir + "/deeptools_frag_lengths.txt",
    qcdir + "/deeptools_frag_lengths.png",
#frag_bamcoverage
    expand(qcdir + "/{library}_bamcoverage.bg",
        library = FRAG_LIBS,
        ),
#frag_plotcoverage
    qcdir + "/frag_coverage.tsv",
    qcdir + "/frag_coverage.pdf",
#frag_multiqc
    qcdir + "/frag_multiqc.html",
    qcdir + "/frag_multiqc_data/multiqc_fastqc.txt",
    qcdir + "/frag_multiqc_data/multiqc_samtools_stats.txt",
    qcdir + "/frag_multiqc_data/multiqc_picard_wgsmetrics.txt",
    qcdir + "/frag_multiqc_data/multiqc_samtools_flagstat.txt",
#frag_make_qc_tsv
    qcdir + "/frag_read_qc.tsv",
    qcdir + "/frag_frag_len.tsv",

rule symlink_inputs:
    input:
        lambda wildcards: frag_lib_dict[wildcards.library],
    output:
        read1 = data_dir + "/analysis/cfdna_wgs/fastqs/{library}_raw_"+ PE1 + ".fastq.gz",
        read2 = data_dir + "/analysis/cfdna_wgs/fastqs/{library}_raw_"+ PE2 + ".fastq.gz",
    params:
        outdir = frag_fastqs,
        script = frag_scriptdir+ "/symlink.sh",
    shell:
        """
        {params.script} \
        {input} \
        {output.read1} \
        {output.read2} \
        {params.outdir}
        """
rule frag_index:
    benchmark: bench_dir + "/frag_index_"+ fasta_base +".benchmark.txt",
    input: genome_fasta,
    log: logdir + "/reference_index/frag_index_"+ fasta_base +".log",
    output: refdir + "/"+ fasta_base +".sa",
    params:
        out_prefix = refdir + "/" + fasta_base,
        script = frag_scriptdir+ "/frag_index.sh",
        log_dir = logdir + "/reference_index",
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} {input} {params.out_prefix} &> {log}
        """
# Adapter-trim and QC reads with fastp
rule frag_fastp:
    benchmark: bench_dir + "/{library}_frag_fastp.benchmark.txt",
    input:
        read1 = frag_fastqs + "/{library}_raw_"+ PE1 + ".fastq.gz",
        read2 = frag_fastqs + "/{library}_raw_"+ PE2 + ".fastq.gz", 
    log:
        cmd = logdir + "/fastp/{library}_frag_fastp.log",
        html = logdir + "/fastp/{library}_frag_fastp.html",
        json = logdir + "/fastp/{library}_frag_fastp.json",
    output:
        read1 = frag_fastqs + "/{library}_processed_"+ PE1 + ".fastq.gz",
        read2 = frag_fastqs + "/{library}_processed_"+ PE2 + ".fastq.gz",
        failed = frag_fastqs + "/{library}_failed_fastp.fastq.gz",
        unpaired1 = frag_fastqs + "/{library}_unpaired_"+ PE1 + ".fastq.gz",
        unpaired2 = frag_fastqs + "/{library}_unpaired_"+ PE2 + ".fastq.gz",
    params:
        script = frag_scriptdir+ "/fastp.sh",
        threads = frag_threads,
        tmp = config.get("tmp", "/tmp"),
        log_dir = logdir + "/fastp",
    resources:
        mem_mb = 180000,
        tmpdir = config.get("tmp", "/tmp"),
    shell:
        """
        mkdir -p {params.log_dir}
        module load fastp
        fastp --detect_adapter_for_pe \
            --failed_out {output.failed} \
            --in1 {input.read1} \
            --in2 {input.read2} \
            --html {log.html} \
            --json {log.json} \
            --out1 {output.read1} \
            --out2 {output.read2} \
            --unpaired1 {output.unpaired1} \
            --unpaired2 {output.unpaired2} \
            --thread {params.threads}\
        > {log.cmd} 2>&1
        """
# Align reads with BWA
rule frag_align:
    benchmark: bench_dir + "/{library}_frag_align.benchmark.txt",
    input:
        ref = genome_fasta,
        read1 = frag_fastqs + "/{library}_processed_"+ PE1 + ".fastq.gz",
        read2 = frag_fastqs + "/{library}_processed_"+ PE2 + ".fastq.gz",
    log: logdir + "/align/{library}_frag_align.log",
    output:
        sort = frag_bams + "/{library}_raw.bam",
        index = frag_bams + "/{library}_raw.bam.bai",
    params:
        script = frag_scriptdir+ "/align.sh",
        threads = config["threads"],
        log_dir = logdir + "/align",
    resources:
        mem_mb = 180000,
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input.ref} \
        {input.read1} \
        {input.read2} \
        {params.threads} \
        {output.sort} &> {log}
        """
# Remove PCR duplicates from aligned reads
rule frag_dedup:
    benchmark: bench_dir + "/{library}_frag_dedup.benchmark.txt",
    input: frag_bams + "/{library}_raw.bam",
    log: logdir + "/dedup/{library}_frag_dedup.log",
    output: frag_bams + "/{library}_dedup.bam",
    params:
        script = frag_scriptdir + "/dedup.sh",
        threads = frag_threads,
        log_dir = logdir + "/dedup",
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {output} \
        {params.threads} &> {log}
        """
# Filter de-duplicated alignments.
# Remove unmapped, not primary, and duplicate reads. Additional location filter by config bedfile variable.
checkpoint frag_filter_alignment:
    benchmark: bench_dir + "/{library}_frag_filter_alignment.benchmark.txt",
    input: frag_bams + "/{library}_dedup.bam",
    log: logdir + "/filter/{library}_frag_filter_alignment.log",
    output: frag_bams + "/{library}_filt.bam",
    params:
        script = frag_scriptdir+ "/filter_alignment.sh",
        threads = frag_threads,
        log_dir = logdir + "/filter",
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """
# Get read quality by FASTQC
rule frag_fastqc:
    benchmark: bench_dir + "/{library}_{processing}_{read}_frag_fastqc.benchmark.txt",
    input: frag_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log: logdir + "/fastqc/{library}_{processing}_{read}_frag_fastqc.log",
    output:
        qcdir + "/{library}_{processing}_{read}_fastqc.html",
        qcdir + "/{library}_{processing}_{read}_fastqc.zip",
    params:
        outdir = qcdir,
        script = frag_scriptdir+ "/fastqc.sh",
        threads = frag_threads,
        log_dir = logdir + "/fastqc"
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """
# Get alignment QC using samtools
rule frag_alignment_qc:
    input: frag_bams + "/{library}_{processing}.bam",
    log:
        flagstat = logdir + "/{library}_{processing}_flagstat_frag_alignment_qc.log",
        samstat = logdir + "/{library}_{processing}_samstats_frag_alignment_qc.log",
    output:
        flagstat = qcdir + "/{library}_{processing}_flagstat.txt",
        samstat = qcdir + "/{library}_{processing}_samstats.txt",
    params:
        script = frag_scriptdir+ "/alignment_qc.sh",
        threads = frag_threads,
    resources:
        time="96:00:00",
    shell:
        """
        {params.script} \
        {input} \
        {log.flagstat} \
        {log.samstat} \
        {output.flagstat} \
        {output.samstat} \
        {params.threads} &> {log}
        """
# Sequencing depth metrics via Picard
rule frag_picard_depth:
    benchmark: bench_dir + "/{library}_frag_picard_depth.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_frag_picard_depth.log",
    output: qcdir + "/{library}_picard_depth.txt",
    params:
        script = frag_scriptdir+ "/picard_depth.sh",
        threads = frag_threads,
    resources:
        mem_mb = 180000,
        cpus = 40,
    shell:
        """
        {params.script} \
        {input} \
        {config[genome_fasta]} \
        {output} \
        {params.threads} &> {log}
        """
# Get fragment sizes using deepTools
rule frag_bampefragsize:
    benchmark: bench_dir + "/frag_bampefragsize.benchmark.txt",
    input: expand(frag_bams + "/{library}_filt.bam", library = FRAG_LIBS),
    log: logdir + "/frag_bampefragsize.log",
    output:
        raw = qcdir + "/deeptools_frag_lengths.txt",
        hist = qcdir + "/deeptools_frag_lengths.png",
    params:
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/bampefragsize.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {log} \
        {output.hist} \
        {output.raw} \
        {params.blacklist} \
        {params.threads} &> {log}
        """
# Make deeptools bamCoverage bedfile
rule frag_bamcoverage:
    benchmark: bench_dir + "/{library}_frag_bamcoverage.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_frag_bamcoverage.log",
    output: qcdir + "/{library}_bamcoverage.bg",
    params:
        bin = "10000",
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/bamcoverage.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.bin} \
        {params.blacklist} \
        {params.threads} &> {log}
        """
# Make deepTools plotCoverage coverage maps for all filtered bams
rule frag_plotcoverage:
    benchmark: bench_dir + "/frag_plotcoverage.benchmark.txt",
    resources:
        time="3:00:00"
    input: expand(frag_bams + "/{library}_filt.bam", library = FRAG_LIBS),
    log: logdir + "/frag_plotcoverage.log",
    output:
        raw = qcdir + "/frag_coverage.tsv",
        plot = qcdir + "/frag_coverage.pdf",
    params:
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/plotcoverage.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.blacklist} \
        {params.threads} \
        {output.raw} \
        {output.plot} &> {log}
        """
# Aggregate QC files using MultiQC
rule frag_multiqc:
    benchmark: bench_dir + "/frag_multiqc.benchmark.txt",
    resources:
        time="3:00:00"
    input:
        expand(logdir + "/fastp/{library}_frag_fastp.json", library = FRAG_LIBS),
        expand(qcdir + "/{library}_{processing}_{read}_fastqc.zip", library = FRAG_LIBS, processing = ["raw", "processed", "unpaired"], read = [PE1,PE2]),
        expand(qcdir + "/{library}_{processing}_samstats.txt", library = FRAG_LIBS, processing = ["raw","filt"]),
        expand(qcdir + "/{library}_{processing}_flagstat.txt", library = FRAG_LIBS, processing = ["raw","filt"]),
        expand(qcdir + "/{library}_picard_depth.txt", library = FRAG_LIBS),
        qcdir + "/deeptools_frag_lengths.txt",
        qcdir + "/frag_coverage.tsv",
    log: logdir + "/frag_multiqc.log",
    output:
        qcdir + "/frag_multiqc.html",
        qcdir + "/frag_multiqc_data/multiqc_fastqc.txt",
        qcdir + "/frag_multiqc_data/multiqc_samtools_stats.txt",
        qcdir + "/frag_multiqc_data/multiqc_picard_wgsmetrics.txt",
        qcdir + "/frag_multiqc_data/multiqc_samtools_flagstat.txt",
    params:
        out_dir = qcdir,
        out_name = "frag_multiqc",
        script = frag_scriptdir+ "/multiqc.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.out_name} \
        {params.out_dir} &> {log}
        """
# Make a tab-separated aggregate QC table
checkpoint frag_make_qc_tsv:
    benchmark: bench_dir + "/frag_make_qc_tsv.benchmark.txt",
    input:
        fq = qcdir + "/frag_multiqc_data/multiqc_fastqc.txt",
        mqsam = qcdir + "/frag_multiqc_data/multiqc_samtools_stats.txt",
        mqflag = qcdir + "/frag_multiqc_data/multiqc_samtools_flagstat.txt",
        picard = qcdir + "/frag_multiqc_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = qcdir + "/deeptools_frag_lengths.txt",
        deeptools_cov = qcdir + "/frag_coverage.tsv",
    log: logdir + "/frag_make_qc_tsv.log",
    output:
        readqc = qcdir + "/frag_read_qc.tsv",
        fraglen = qcdir + "/frag_frag_len.tsv",
    params:
        script = frag_scriptdir+ "/make_qc_tsv.R",
    shell:
        """
        module load R
        Rscript {params.script} \
        {input.fq} \
        {input.mqsam} \
        {input.mqflag} \
        {input.picard} \
        {input.deeptools_frag} \
        {input.deeptools_cov} \
        {output.readqc} \
        {output.fraglen} >& {log}
        """

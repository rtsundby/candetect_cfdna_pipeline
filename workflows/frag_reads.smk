##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np
import os

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
bench_dir           = data_dir + "/benchmark"
frag_analysis_dir   = data_dir + "/analysis/cfdna_wgs"
frag_bams           = data_dir  + "/bams"
frag_fastqs         = data_dir  + "/fastq"
frag_scriptdir      = frag_repo +  "/scripts"
logdir              = data_dir + "/logs/preprocessing"
qcdir               = data_dir  + "/qc"

#Make sure logdir exists
os.makedirs(logdir, exist_ok=True)

####------ Read samples ------####
frag_libs = pd.read_table(config["samples"])

# Add detailed logging about frag_libs
print("\n=== frag_libs DataFrame Info ===")
print("Number of libraries:", len(frag_libs))
print("\nLibrary names:", frag_libs["library"].tolist())
print("\nCohorts:", frag_libs["cohort"].unique().tolist())
print("\nFile paths:", frag_libs["file"].tolist())
print("\nLibrary types:", frag_libs["library_type"].unique().tolist())
print("\nIsolation types:", frag_libs["isolation_type"].unique().tolist())
print("\n=== End frag_libs Info ===\n")

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

# Debug print
print("Available library names:", FRAG_LIBS)
print("Library dictionary keys:", list(frag_lib_dict.keys()))

# Make  a list of healthy libraries
FRAG_HEALTHY_LIBRARIES = frag_libs[frag_libs['cohort'] == normal_control]['library'].tolist()

rule all:
  input: 
#symlinks
    expand(frag_fastqs + "/{library}_raw_"+ PE1 + ".fastq.gz",
            library = FRAG_LIBS,
    ),
    expand(frag_fastqs + "/{library}_raw_"+ PE2 + ".fastq.gz",
            library = FRAG_LIBS,
    ),
#frag_index
    refdir + "/"+ fasta_base +".sa",
#frag_fastp
    expand(frag_fastqs + "/{library}_processed_"+ PE1 + ".fastq.gz",
        library = FRAG_LIBS,
    ),
    expand(frag_fastqs + "/{library}_processed_"+ PE2 + ".fastq.gz",
       library = FRAG_LIBS,
    ),
#frag_align
    expand(frag_bams + "/{library}_raw.bam",
        library = FRAG_LIBS,
    ),
    expand(frag_bams + "/{library}_raw.bam.bai",
        library = FRAG_LIBS,
    ),
#frag_dedup
    expand(frag_bams + "/{library}_dedup.bam",
        library = FRAG_LIBS,
    ),
#frag_filter
    expand(frag_bams + "/{library}_filt.bam",
        library = FRAG_LIBS,
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
    expand(qcdir + "/{library}_coverage.tsv",
        library = FRAG_LIBS,
        ),
    expand(qcdir + "/{library}_coverage.pdf",
        library = FRAG_LIBS,
        ),
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
    resources:
        time   = "00:10:00",
        mem_gb = "16g",
        cpus   = "2",
    input:
        lambda wildcards: frag_libs[frag_libs["library"] == wildcards.library]["file"],
    output:
        read1 = frag_fastqs + "/{library}_raw_"+ PE1 + ".fastq.gz",
        read2 = frag_fastqs + "/{library}_raw_"+ PE2 + ".fastq.gz",
    params:
        outdir = frag_fastqs,
        script = frag_scriptdir+ "/detect_read_pairs.sh",
        log_dir = logdir + "/symlink",
        pair_detection_dir = logdir + "/pair_detection",
    log:
        cmd = logdir + "/symlink/{library}_symlink.log",
        debug = logdir + "/symlink/{library}_symlink_debug.log",
        pair_detection = logdir + "/pair_detection/{library}_read_pair_detection.log",
    run:
        # Create log directories
        os.makedirs(params.log_dir, exist_ok=True)
        os.makedirs(params.pair_detection_dir, exist_ok=True)
        
        # Get source file path
        source_file = frag_libs[frag_libs["library"] == wildcards.library]["file"].iloc[0]
        
        # Log debug information
        with open(log.debug, 'w') as debug_log:
            debug_log.write(f"=== Symlink Debug Info ===\n")
            debug_log.write(f"Library: {wildcards.library}\n")
            debug_log.write(f"Source file: {source_file}\n")
            debug_log.write(f"Output PE1: {output.read1}\n")
            debug_log.write(f"Output PE2: {output.read2}\n")
            debug_log.write(f"Output directory: {params.outdir}\n")
            debug_log.write(f"File exists: {os.path.exists(source_file)}\n")
            debug_log.write(f"File readable: {os.access(source_file, os.R_OK)}\n")
            debug_log.write(f"All libraries in FRAG_LIBS: {FRAG_LIBS}\n")
            debug_log.write(f"Current library in wildcards: {wildcards.library}\n")
            debug_log.write(f"Library exists in FRAG_LIBS: {wildcards.library in FRAG_LIBS}\n")
            debug_log.write(f"=== End Debug Info ===\n")
        
        # Check if source file exists
        if not os.path.exists(source_file):
            error_msg = f"Source file not found: {source_file}"
            with open(log.debug, 'a') as debug_log:
                debug_log.write(f"ERROR: {error_msg}\n")
            raise FileNotFoundError(error_msg)
            
        # Create output directory
        os.makedirs(params.outdir, exist_ok=True)
        
        # Create symlinks
        shell("""
            echo "Creating symlinks for {wildcards.library}" > {log.cmd}
            echo "Source: {input}" >> {log.cmd}
            echo "Destination R1: {output.read1}" >> {log.cmd}
            echo "Destination R2: {output.read2}" >> {log.cmd}
            
            {params.script} \
            {input} \
            {output.read1} \
            {output.read2} \
            {params.outdir} \
            >> {log.pair_detection} 2>&1
            
            echo "Symlink creation completed" >> {log.cmd}
            ls -l {output.read1} {output.read2} >> {log.cmd} 2>&1
        """)
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
            --disable_length_filtering \
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
        time = "72:00:00",
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
        flagstat = logdir + "/alignment_qc/{library}_{processing}_flagstat.log",
        samstat = logdir + "/alignment_qc/{library}_{processing}_samstats.log",
    output:
        flagstat = qcdir + "/{library}_{processing}_flagstat.txt",
        samstat = qcdir + "/{library}_{processing}_samstats.txt",
    params:
        script = frag_scriptdir+ "/alignment_qc.sh",
        threads = frag_threads,
        log_dir = logdir + "/alignment_qc",
    resources:
        time="96:00:00",
    shell:
        """
        mkdir -p {params.log_dir}
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
    log: logdir + "/picard_depth/{library}_picard_depth.log",
    output: qcdir + "/{library}_picard_depth.txt",
    params:
        script = frag_scriptdir+ "/picard_depth.sh",
        threads = frag_threads,
        log_dir = logdir + "/picard_depth",
    resources:
        mem_mb = 180000,
        cpus = 40,
    shell:
        """
        mkdir -p {params.log_dir}
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
    log: logdir + "/bampefragsize/bampefragsize.log",
    output:
        raw = qcdir + "/deeptools_frag_lengths.txt",
        hist = qcdir + "/deeptools_frag_lengths.png",
    params:
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/bampefragsize.sh",
        threads = frag_threads,
        log_dir = logdir + "/bampefragsize",
    shell:
        """
        mkdir -p {params.log_dir}
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
    log: logdir + "/bamcoverage/{library}_bamcoverage.log",
    output: qcdir + "/{library}_bamcoverage.bg",
    params:
        bin = "10000",
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/bamcoverage.sh",
        threads = frag_threads,
        read_type = config.get("seq_params", {}).get("read_type", "PE"),
        read_length = config.get("seq_params", {}).get("read_length", 150),
        log_dir = logdir + "/bamcoverage",
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        {input} \
        {output} \
        {params.bin} \
        {params.blacklist} \
        {params.threads} \
        {params.read_type} \
        {params.read_length} &> {log}
        """
# Make deepTools plotCoverage coverage maps for all filtered bams
rule frag_plotcoverage:
    benchmark: bench_dir + "/{library}_frag_plotcoverage.benchmark.txt",
    resources:
        time="3:00:00"
    input: 
        frag_bams + "/{library}_filt.bam"
    log: logdir + "/plotcoverage/{library}_plotcoverage.log",
    output:
        raw = qcdir + "/{library}_coverage.tsv",
        plot = qcdir + "/{library}_coverage.pdf",
    params:
        blacklist = config["blklist"],
        script = frag_scriptdir+ "/plotcoverage.sh",
        threads = frag_threads,
        log_dir = logdir + "/plotcoverage",
    shell:
        """
        mkdir -p {params.log_dir}
        {params.script} \
        "{input}" \
        {params.blacklist} \
        {params.threads} \
        {output.raw} \
        {output.plot} &> {log}
        """
# Aggregate QC files using MultiQC
rule frag_multiqc:
    benchmark: bench_dir + "/frag_multiqc.benchmark.txt"
    resources:
        time="3:00:00"
    input:
        fastp = expand(logdir + "/fastp/{library}_frag_fastp.json", library=FRAG_LIBS),  
        fastqc = expand(qcdir + "/{library}_{processing}_{read}_fastqc.zip", 
                       library=FRAG_LIBS, 
                       processing=["raw", "processed", "unpaired"], 
                       read=[PE1,PE2]),
        samstats = expand(qcdir + "/{library}_{processing}_samstats.txt", 
                         library=FRAG_LIBS, 
                         processing=["raw","filt"]),
        flagstat = expand(qcdir + "/{library}_{processing}_flagstat.txt", 
                         library=FRAG_LIBS, 
                         processing=["raw","filt"]),
        picard = expand(qcdir + "/{library}_picard_depth.txt", library=FRAG_LIBS),
        frag_len = qcdir + "/deeptools_frag_lengths.txt",
        coverage = expand(qcdir + "/{library}_coverage.tsv", library=FRAG_LIBS)
    output:
        html = qcdir + "/frag_multiqc.html",
        fastqc = qcdir + "/frag_multiqc_data/multiqc_fastqc.txt",
        samstats = qcdir + "/frag_multiqc_data/multiqc_samtools_stats.txt",
        picard = qcdir + "/frag_multiqc_data/multiqc_picard_wgsmetrics.txt",
        flagstat = qcdir + "/frag_multiqc_data/multiqc_samtools_flagstat.txt"
    log: logdir + "/multiqc/multiqc.log"
    params:
        script = frag_scriptdir + "/multiqc.sh",
        input_files = lambda w, input: " ".join([
            *input.fastp,
            *input.fastqc,
            *input.samstats,
            *input.flagstat,
            *input.picard,
            input.frag_len,
            *input.coverage
        ]),
        log_dir = logdir + "/multiqc",
    shell:
        """
        mkdir -p {qcdir}/frag_multiqc_data
        mkdir -p {params.log_dir}
        
        {params.script} \
            '{params.input_files}' \
            frag_multiqc \
            {qcdir} \
            &> {log}
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
        deeptools_cov = expand(qcdir + "/{library}_coverage.tsv", library = FRAG_LIBS)
    log: logdir + "/make_qc_tsv/make_qc_tsv.log",
    output:
        readqc = qcdir + "/frag_read_qc.tsv",
        fraglen = qcdir + "/frag_frag_len.tsv",
        deeptools_cov = qcdir +"/frag_coverage.tsv",
    params:
        script = frag_scriptdir+ "/make_qc_tsv.R",
        script_cov = frag_scriptdir + "/combine_deeptool_coverage.R",
        log_dir = logdir + "/make_qc_tsv",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {qcdir}/summary
        module load R

        Rscript {params.script_cov} --out {output.deeptools_cov} {input.deeptools_cov}

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

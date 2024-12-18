#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

import pandas as pd
import re
import numpy as np

# Values directly from configuration file
samples = config["samples"]
metadata_file= config["metadata_file"]
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

# Directory values derived from data_dir in configuration YAML
data_dir             = config["data_dir"]
frag                 = data_dir + "/analysis/cfdna_wgs/frag"
frag_bams            = config["frag_bams"]
frag_fastqs          = config["raw_fastq"]
frag_frag            = data_dir + "/analysis/cfdna_wgs/frag/frag"
frag_frag_beds       = data_dir + "/analysis/cfdna_wgs/frag/frag/beds"
frag_frag_counts     = data_dir + "/analysis/cfdna_wgs/frag/frag/counts"
frag_frag_gc_distros = data_dir + "/analysis/cfdna_wgs/frag/frag/distros"
qcdir                = data_dir + "/analysis/cfdna_wgs/qc"
benchdir             = data_dir + "/benchmark"
logdir               = data_dir + "/logs"
refdir               = config["ref_dir"]
bwa_dir              = config["ref_dir"]
motif_dir            = data_dir + "/analysis/cfdna_wgs/motifs"
frag_scriptdir = config["frag_repo"] +  "/scripts"

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
#bam_to_motif
    expand(motif_dir + "/{library}_motifs.tsv",
        library = frag_libs["library"]),
#motif_merge
    data_dir + "/analysis/cfdna_wgs/features/combined_motif_counts_matrix.tsv",
    data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
#mds
    data_dir + "/analysis/cfdna_wgs/features/mds.RData",
    data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv",
#fprofile
    data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
    data_dir + "/analysis/cfdna_wgs/features/motif_per_fprofile.tsv",
#fprofile plot
    data_dir + "/analysis/cfdna_wgs/features/fprofiles.pdf",
#fprofile_roc_plots:
    data_dir + "/analysis/cfdna_wgs/features/fprofile_ROC.pdf",
    data_dir + "/analysis/cfdna_wgs/features/fprofile_AUC.csv",
rule bam_to_motif:
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/motifs/bam_to_motif/{library}_bam_to_motif.log",
    output: motif_dir + "/{library}_motifs.tsv",
    params:
        fasta = genome_fasta,
        script = frag_scriptdir + "/sample_motifs.py",
        threads = frag_threads,
        out_dir = motif_dir,
        log_dir = logdir + "/motifs/bam_to_motif",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        module load python samtools bedtools
        {params.script} --bam_file {input} --output_file {output} --motif_length 4 \
        --threads {params.threads} --num_reads 100000 --reference_genome {params.fasta} \
        > {log} 2>&1
        """
# Merge short and long fragment counts by genomic poistion for all libraries
rule motif_merge:
    input: expand(motif_dir + "/{library}_motifs.tsv",  library = FRAG_LIBS),
    log: logdir + "/motifs/motif_merge/motif_merge.log",
    output:  
        combined_counts = data_dir + "/analysis/cfdna_wgs/features/combined_motif_counts_matrix.tsv",
        combined_fractions = data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    params:
        motif_tsv_dir = motif_dir,
        script = frag_scriptdir + "/motif_data_count.R",
        out_dir = data_dir + "/analysis/cfdna_wgs/features",
        log_dir = logdir + "/motifs/motif_merge",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        module load R
        Rscript {params.script} \
        {params.motif_tsv_dir} \
        {output.combined_counts} \
        {output.combined_fractions} \
        &> {log}
        """
rule mds:
    input: data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    log: logdir + "/motifs/mds/motif_diversity_score.log",
    output:  
        mdsRDS = data_dir + "/analysis/cfdna_wgs/features/mds.RData",
        CSVout = data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv",
    params:
        annotation_table = config["samples"],
        script = frag_scriptdir + "/motif_diversity_score.R",
        out_dir = data_dir + "/analysis/cfdna_wgs/features",
        log_dir = logdir + "/motifs/mds",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        module load R
        Rscript {params.script} \
        {input} \
        {params.annotation_table} \
        {output.mdsRDS} \
        {output.CSVout} \
        &> {log}
        """
rule f_profile:
    input: data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    log: logdir + "/motifs/f_profile/f-profile.log",
    output: 
        fprofile_per_lib = data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
        motif_per_fprofile = data_dir + "/analysis/cfdna_wgs/features/motif_per_fprofile.tsv",
        model_dir = directory(data_dir + "/analysis/cfdna_wgs/features/fprofile_model/"),
    params:
        script = frag_scriptdir + "/fprofiles_pkl.py",
        log_dir = logdir + "/motifs/f_profile",
        existing_model = config["existing_motif_model"],
    run:
        import os
        import shutil
        os.makedirs(output.model_dir, exist_ok=True)

        # If loading an existing model, create symlink for pkl into output model_dir
        if params.existing_model != "no":
            model_file = os.path.join(params.existing_model, "fprofile_model.pkl")
            symlink_target = os.path.join(output.model_dir, "fprofile_model.pkl")

            # Create a symlink to the existing model if it doesn't already exist
            if not os.path.exists(symlink_target):
                if os.path.exists(model_file):
                    shutil.copy(model_file, symlink_target)
                    print(f"Symlink created for existing model: {symlink_target}")
                else:
                    raise FileNotFoundError(f"Specified model file not found: {model_file}")

            # Use the existing model for analysis
            shell(f"""
            mkdir -p {params.log_dir}
            module load python
            {params.script} --data_file {input.data_file} \
            --output_fprof_per_lib {output.fprofile_per_lib} \
            --output_motif_per_fprof {output.motif_per_fprofile} \
            --load_model {symlink_target} \
            > {log} 2>&1
            """)
        else:
            # Train a new model and save it in the output model_dir
            model_output = os.path.join(output.model_dir, "fprofile_model.pkl")
            shell(f"""
            mkdir -p {params.log_dir}
            module load python
            {params.script} --data_file {input.data_file} \
            --output_fprof_per_lib {output.fprofile_per_lib} \
            --output_motif_per_fprof {output.motif_per_fprofile} \
            --save_model {model_output} \
            > {log} 2>&1
            """)
#    shell:
#        """
#        mkdir -p {params.log_dir}
#        module load python
#        {params.script} --data_file {input} \
#        --output_fprof_per_lib {output.fprofile_per_lib} \
#        --output_motif_per_fprof {output.motif_per_fprofile} \
#        > {log} 2>&1
#        """
rule f_profile_plot:
    input: data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
    log: logdir + "/f-profile_plot.log",
    output: 
        data_dir + "/analysis/cfdna_wgs/features/fprofiles.pdf",
    params:
        script = frag_scriptdir + "/fprofile_stackedbarplot.py",
        plot_title = config["fprofile_stacked_title"],
        metadata_file = metadata_file, #Note, if not using e.g. different cycles, set as "none" in config
        colorscheme = "viridis",
    shell:
        """
        module load python
        python {params.script} --input_file {input} --metadata_file {params.metadata_file} \
        --output_file {output} --title "{params.plot_title}" --colormap {params.colorscheme} \
        > {log} 2>&1
	    """      
rule fprofile_roc_plots:
    input: data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
    log: logdir + "motifs/roc_plots/roc_f-profiles.log",
    output:  
        plot_out = data_dir + "/analysis/cfdna_wgs/features/fprofile_ROC.pdf",
        csv_out = data_dir + "/analysis/cfdna_wgs/features/fprofile_AUC.csv",
    params:
        libraries_tsv = config["samples"],
        script = frag_scriptdir + "/motif_ROC_plot.R",
        out_dir = data_dir + "/analysis/cfdna_wgs/features",
        log_dir = logdir + "/motifs/roc_plots",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        module load R
        Rscript {params.script} \
        {input} \
        {params.libraries_tsv} \
        {output.plot_out} \
        {output.csv_out} \
        &> {log}
        """
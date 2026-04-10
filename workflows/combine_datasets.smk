#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Create a master dataframe with fragmentomic outputs (R version)          #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

import os
import itertools as it

# Define directories
data_dir = config["data_dir"]
log_dir = data_dir + "/logs/combine_datasets"
plot_dir = data_dir + "/analysis/cfdna_wgs/plots"

# Get cohorts and normal control from config
COHORTS = config["cohorts"]
NORMAL_CONTROL = config["normal_control"]

# Create directories if they don't exist
for d in [log_dir, plot_dir]:
    os.makedirs(d, exist_ok=True)

# Define cohort comparisons (one-versus-one)
COHORT_1 = []
COHORT_2 = []
for c1, c2 in it.combinations(COHORTS,2):
    COHORT_1.append(c1)
    COHORT_2.append(c2)
    print(c1,c2)

rule all:
    input:
        expand(data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv",
            zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        expand(plot_dir + "/subject_trajectories_{c1}_v_{c2}.png",
            zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        expand(plot_dir + "/auc_comparison_{c1}_v_{c2}.png",
            zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        ),
        expand(plot_dir + "/distributions_{c1}_v_{c2}",
            zip,
            c1 = COHORT_1,
            c2 = COHORT_2
        )

rule combine_datasets:
    resources:
        mem_mb = "64000",
        time   = "0:30:00",
        cpus   = "16"
    input:
        tf = lambda wildcards: data_dir + "/analysis/cfdna_wgs/ichorCNA/tf.csv" if os.path.exists(data_dir + "/analysis/cfdna_wgs/ichorCNA/tf.csv") else [],
        fragle = lambda wildcards: data_dir + "/FRAGLE/master_compiled_csv_2025-04-22.csv" if os.path.exists(data_dir + "/FRAGLE/master_compiled_csv_2025-04-22.csv") else [],
        #fragle = lambda wildcards: data_dir + "/FRAGLE/master_compiled_csv_2025-04-13.csv" if os.path.exists(data_dir + "/FRAGLE/master_compiled_csv_2025-04-13.csv") else [],
        mds = lambda wildcards: data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv" if os.path.exists(data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv") else [],
        fprofile_auc = lambda wildcards: data_dir + "/analysis/cfdna_wgs/features/fprofile_AUC.csv" if os.path.exists(data_dir + "/analysis/cfdna_wgs/features/fprofile_AUC.csv") else [],
        fprofiles = lambda wildcards: data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv" if os.path.exists(data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv") else [],
        nmf = lambda wildcards: data_dir + "/analysis/cfdna_wgs/nmf/aucs/aucs.csv" if os.path.exists(data_dir + "/analysis/cfdna_wgs/nmf/aucs/aucs.csv") else [],
        binwise = data_dir + "/analysis/cfdna_wgs/bin_wise_fragmentomics/model/probs/{c1}_v_{c2}.csv"
    output:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv",
        status_file = data_dir + "/analysis/cfdna_wgs/combined_dataset_status_{c1}_v_{c2}.txt"
    log:
        log_dir + "/combine_datasets_{c1}_v_{c2}.log"
    params:
        tf_arg = lambda w, input: f"--tf='{input.tf}'" if input.tf else "",
        fragle_arg = lambda w, input: f"--fragle='{input.fragle}'" if input.fragle else "",
        mds_arg = lambda w, input: f"--mds='{input.mds}'" if input.mds else "",
        fprofile_auc_arg = lambda w, input: f"--fprofile_auc='{input.fprofile_auc}'" if input.fprofile_auc else "",
        fprofiles_arg = lambda w, input: f"--fprofiles='{input.fprofiles}'" if input.fprofiles else "",
        nmf_arg = lambda w, input: f"--nmf='{input.nmf}'" if input.nmf else "",
        binwise_arg = lambda w, input: f"--binwise='{input.binwise}'" if input.binwise else ""
    shell:
        """
        module load R;
        Rscript scripts/combine_datasets.R \
            --config='{config[samples]}' \
            --data_dir='{data_dir}' \
            --cohort1='{wildcards.c1}' \
            --cohort2='{wildcards.c2}' \
            --normal_control='{NORMAL_CONTROL}' \
            --output='{output.combined_dataset}' \
            --status='{output.status_file}' \
            --log='{log}' \
            {params.tf_arg} \
            {params.fragle_arg} \
            {params.mds_arg} \
            {params.fprofile_auc_arg} \
            {params.fprofiles_arg} \
            {params.nmf_arg} \
            {params.binwise_arg} 2>&1 | tee -a {log}
        """

rule plot_subject_trajectories:
    input:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv"
    output:
        plot = plot_dir + "/subject_trajectories_{c1}_v_{c2}.png"
    log:
        log_dir + "/plot_subject_trajectories_{c1}_v_{c2}.log"
    resources:
        mem_mb = 16000,
        time = "0:10:00",
        cpus = 4,
        lscratch = 100
    shell:
        """
        mkdir -p {plot_dir}
        module load R;
        Rscript scripts/plot_subject_trajectories.R \
            --input {input.combined_dataset} \
            --output {output.plot} \
            --cohort1 {wildcards.c1} \
            --cohort2 {wildcards.c2} \
            --normal_control {NORMAL_CONTROL} \
            --log {log} 2>&1 | tee -a {log}
        """

rule plot_auc_comparison:
    input:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv"
    output:
        plot = plot_dir + "/auc_comparison_{c1}_v_{c2}.png"
    log:
        log_dir + "/plot_auc_comparison_{c1}_v_{c2}.log"
    resources:
        mem_mb = 16000,
        time = "0:10:00",
        cpus = 4,
        lscratch = 100
    shell:
        """
        mkdir -p {plot_dir}
        module load R;
        Rscript scripts/plot_auc.R \
            --input {input.combined_dataset} \
            --output {output.plot} \
            --cohort1 {wildcards.c1} \
            --cohort2 {wildcards.c2} \
            --normal_control {NORMAL_CONTROL} \
            --log {log} 2>&1 | tee -a {log}
        """

rule plot_metric_distributions:
    input:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv"
    output:
        directory(plot_dir + "/distributions_{c1}_v_{c2}")
    log:
        log_dir + "/plot_metric_distributions_{c1}_v_{c2}.log"
    resources:
        mem_mb = 16000,
        time = "0:10:00",
        cpus = 4,
        lscratch = 100
    shell:
        """
        mkdir -p {output}
        module load R;
        Rscript scripts/plot_metric_distributions.R \
            --input {input.combined_dataset} \
            --output_dir {output} \
            --cohort1 {wildcards.c1} \
            --cohort2 {wildcards.c2} \
            --normal_control {NORMAL_CONTROL} \
            --log {log} 2>&1 | tee -a {log}
        """

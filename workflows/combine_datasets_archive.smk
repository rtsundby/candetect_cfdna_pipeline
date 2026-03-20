#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Create a master dataframe with fragmentomic outputs                      #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

import pandas as pd
import re
import numpy as np
import os
import glob
import itertools as it
import sys
print(sys.executable)

# Define directories
data_dir = config["data_dir"]
log_dir = data_dir + "/logs/combine_datasets"
plot_dir = data_dir + "/analysis/cfdna_wgs/plots"

# Get cohorts and normal control from config
COHORTS = config["cohorts"]
NORMAL_CONTROL = config["normal_control"]

# Create directories if they don't exist
os.makedirs(log_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

# Get sample list from config
samples_df = pd.read_csv(config["samples"], sep="\t")
samples_df["library"] = samples_df["library"].str.upper()

# Initialize combined dataframe with library and cohort
combined_df = samples_df[["library", "cohort"]].copy()

# Function to read and merge data
def read_and_merge_data(file_path, id_col, value_cols, id_col_combined="library"):
    try:
        df = pd.read_csv(file_path)
        # Clean up ID column if needed
        if id_col == "Sample_ID":
            df[id_col] = df[id_col].str.strip('"').str.upper()
        elif id_col == "library":
            df[id_col] = df[id_col].str.replace("_motifs", "").str.upper()
        elif id_col == "sample_id":
            df[id_col] = df[id_col].str.upper()
        return df[[id_col] + value_cols]
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return pd.DataFrame(columns=[id_col] + value_cols)

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
        )

rule combine_datasets:
    resources:
        mem_mb = "64000",
        time   = "0:30:00",
        cpus   = "16"
    input:
        # Make inputs optional by using a function to check existence
        _ = lambda wildcards: [f for f in [
            data_dir + "/analysis/cfdna_wgs/ichorCNA/tf.csv",
            data_dir + "/FRAGLE/master_compiled_csv_2025-04-13.csv",
            data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv",
            data_dir + "/analysis/cfdna_wgs/features/fprofile_AUC.csv",
            data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
            data_dir + "/analysis/cfdna_wgs/nmf/aucs/aucs.csv",
            data_dir + "/analysis/cfdna_wgs/bin_wise_frag/model/probs/{c1}_v_{c2}.csv"
        ] if os.path.exists(f)]
    output:
        combined_dataset = data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv",
        status_file = data_dir + "/analysis/cfdna_wgs/combined_dataset_status_{c1}_v_{c2}.txt"
    log:
        log_dir + "/combine_datasets_{c1}_v_{c2}.log"
    run:
        try:
            # Basic setup and debugging
            print("Starting combine_datasets rule", file=open(log[0], 'a'))
            print(f"Wildcards: c1={wildcards.c1}, c2={wildcards.c2}", file=open(log[0], 'a'))
            print(f"Data directory: {data_dir}", file=open(log[0], 'a'))
            
            # Create output directory if it doesn't exist
            os.makedirs(os.path.dirname(output.combined_dataset), exist_ok=True)
            print("Created output directory", file=open(log[0], 'a'))
            
            # Debug: Print available input files
            print("Available input files:", file=open(log[0], 'a'))
            for i, f in enumerate(input._()):
                print(f"Input {i}: {f}", file=open(log[0], 'a'))
            
            # Initialize status tracking
            status = {
                'tf': False,
                'fragle': False,
                'mds': False,
                'fprofile': False,
                'nmf': False,
                'binwise': False
            }
            print("Initialized status tracking", file=open(log[0], 'a'))

            # Get sample list from config
            print(f"Reading samples from config: {config['samples']}", file=open(log[0], 'a'))
            samples_df = pd.read_csv(config["samples"], sep="\t")
            samples_df["library"] = samples_df["library"].str.upper()
            print(f"Read {len(samples_df)} samples", file=open(log[0], 'a'))

            # Initialize combined dataframe with library and cohort
            combined_df = samples_df[["library", "cohort"]].copy()
            print("Initialized combined dataframe", file=open(log[0], 'a'))
            
            # Debug: Print initial dataframe
            print("Initial dataframe shape:", combined_df.shape, file=open(log[0], 'a'))
            print("Initial columns:", combined_df.columns.tolist(), file=open(log[0], 'a'))

            # 1. Merge Tumor Fraction data if available
            if len(input._()) > 0 and os.path.exists(input._()[0]):
                print("Processing Tumor Fraction data...", file=open(log[0], 'a'))
                tf_df = read_and_merge_data(
                    input._()[0],
                    "Sample_ID",
                    ["Tumor_Fraction"]
                )
                combined_df = pd.merge(combined_df, tf_df, left_on="library", right_on="Sample_ID", how="left")
                combined_df = combined_df.drop("Sample_ID", axis=1)
                status['tf'] = True
            else:
                combined_df['Tumor_Fraction'] = 'missing'
                print("Warning: Tumor Fraction file not found", file=open(log[0], 'a'))

            # 2. Merge FRAGLE data if available
            if len(input._()) > 1 and os.path.exists(input._()[1]):
                print("Processing FRAGLE data...", file=open(log[0], 'a'))
                fragle_df = read_and_merge_data(
                    input._()[1],
                    "sample_id",
                    ["lt_ctdna_burden"]
                )
                combined_df = pd.merge(combined_df, fragle_df, left_on="library", right_on="sample_id", how="left")
                combined_df = combined_df.drop("sample_id", axis=1)
                status['fragle'] = True
            else:
                combined_df['lt_ctdna_burden'] = 'missing'
                print("Warning: FRAGLE file not found", file=open(log[0], 'a'))

            # 3. Merge motif diversity scores if available
            if input._()[2]:
                mds_df = read_and_merge_data(
                    input._()[2],
                    "library",
                    ["mds"]
                )
                combined_df = pd.merge(combined_df, mds_df, on="library", how="left")
                status['mds'] = True
            else:
                combined_df['mds'] = 'missing'
                print("Warning: Motif diversity scores file not found")

            # 4. Get best f-profiles and merge if available
            if input._()[3] and input._()[4]:
                fprofile_auc = pd.read_csv(input._()[3])
                top_3_fprofiles = fprofile_auc.nlargest(3, "AUC")["fprofile"].tolist()

                fprofiles_df = pd.read_csv(input._()[4], sep="\t")
                fprofiles_df = fprofiles_df[fprofiles_df.index.isin(top_3_fprofiles)]
                fprofiles_df = fprofiles_df.T.reset_index()
                fprofiles_df["index"] = fprofiles_df["index"].str.replace("_motifs", "").str.upper()
                fprofiles_df = fprofiles_df.rename(columns={"index": "library"})

                combined_df = pd.merge(combined_df, fprofiles_df, on="library", how="left")
                status['fprofile'] = True
            else:
                print("Warning: F-profile files not found")
                # Add missing columns for top 3 f-profiles
                for i in range(1, 4):
                    combined_df[f'fprofile_{i}'] = 'missing'

            # 5. Get best NMF components and merge if available
            if input._()[5]:
                aucs_df = pd.read_csv(input._()[5])
                # Filter for current cohort comparison
                cohort_aucs = aucs_df[(aucs_df['cohort1'] == wildcards.c1) & (aucs_df['cohort2'] == wildcards.c2)]
                best_n_components = cohort_aucs.loc[cohort_aucs["med_score"].idxmax(), "n_components"]
                
                print(f"Best number of components for {wildcards.c1} vs {wildcards.c2}: {best_n_components}")
                
                # Find the correct NMF scores file
                nmf_scores_pattern = os.path.join(data_dir + "/analysis/cfdna_wgs/nmf/scores", f"{wildcards.c1}_vs_{wildcards.c2}.n={best_n_components}.csv")
                nmf_scores_files = glob.glob(nmf_scores_pattern)
                
                if nmf_scores_files:
                    print(f"Using NMF scores file: {nmf_scores_files[0]}")
                    
                    nmf_scores_df = read_and_merge_data(
                        nmf_scores_files[0],
                        "library",
                        ["med_score"]
                    )
                    combined_df = pd.merge(combined_df, nmf_scores_df, on="library", how="left")
                    combined_df = combined_df.rename(columns={"med_score": "nmf"})
                    status['nmf'] = True
                else:
                    print(f"Warning: No NMF scores file found for n_components={best_n_components}")
                    combined_df['nmf'] = 'missing'
            else:
                print("Warning: NMF AUCs file not found")
                combined_df['nmf'] = 'missing'

            # 6. Merge bin-wise fragmentomics data if available
            if input._()[6]:
                binwise_df = pd.read_csv(input._()[6])
                binwise_df = binwise_df.rename(columns={"prob": "binwise_prob"})
                combined_df = pd.merge(combined_df, binwise_df, on="library", how="left")
                status['binwise'] = True
            else:
                print(f"Warning: Bin-wise probabilities file not found for {wildcards.c1} vs {wildcards.c2}")
                combined_df['binwise_prob'] = 'missing'

            # Replace NaN with "missing"
            combined_df = combined_df.fillna("missing")

            # Save status information
            with open(output.status_file, 'w') as f:
                f.write("Dataset combination status:\n")
                for key, value in status.items():
                    f.write(f"{key}: {'Included' if value else 'Missing'}\n")

            # Save combined dataset
            combined_df.to_csv(output.combined_dataset, index=False)

        except Exception as e:
            print(f"Error in combine_datasets rule: {str(e)}", file=open(log[0], 'a'))
            print(f"Error type: {type(e).__name__}", file=open(log[0], 'a'))
            import traceback
            print("Traceback:", file=open(log[0], 'a'))
            print(traceback.format_exc(), file=open(log[0], 'a'))
            raise

rule plot_subject_trajectories:
    resources:
        mem_mb = "16000",
        time = "0:10:00",
        cpus = "4"
    input:
        combined_dataset = rules.combine_datasets.output.combined_dataset
    output:
        plot = plot_dir + "/subject_trajectories_{c1}_v_{c2}.png"
    log:
        log_dir + "/plot_subject_trajectories_{c1}_v_{c2}.log"
    run:
        import seaborn as sns
        import matplotlib.pyplot as plt
        # Read the combined dataset
        df = pd.read_csv(input.combined_dataset)
        
        # Filter for the current cohort comparison, excluding healthy controls if they're part of the comparison
        if NORMAL_CONTROL in [wildcards.c1, wildcards.c2]:
            # If healthy is one of the cohorts, only include the other cohort
            other_cohort = wildcards.c2 if wildcards.c1 == NORMAL_CONTROL else wildcards.c1
            df = df[df['cohort'] == other_cohort]
        else:
            # If neither cohort is healthy, include both
            df = df[df['cohort'].isin([wildcards.c1, wildcards.c2])]
        
        # Extract subject and cycle information
        df['subject'] = df['library'].str.split('_CY').str[0]
        df['cycle'] = df['library'].str.extract('_CY(\d+)').astype(int)
        
        # Get all metrics including F-profiles
        base_metrics = ['Tumor_Fraction', 'lt_ctdna_burden', 'mds', 'nmf', 'binwise_prob']
        fprofile_columns = [col for col in df.columns if col.startswith('fprofile_')]
        all_metrics = base_metrics + fprofile_columns
        
        # Convert metrics to numeric, replacing 'missing' with NaN
        for metric in all_metrics:
            df[metric] = pd.to_numeric(df[metric].replace('missing', pd.NA), errors='coerce')
        
        # Create a color palette with a unique color for each subject
        subjects = df['subject'].unique()
        colors = sns.color_palette('husl', n_colors=len(subjects))
        subject_colors = dict(zip(subjects, colors))
        
        # Calculate number of rows needed (2 metrics per row)
        n_metrics = len(all_metrics)
        n_rows = (n_metrics + 1) // 2
        
        # Create the faceted plot
        fig, axes = plt.subplots(n_rows, 2, figsize=(15, 5*n_rows))
        fig.suptitle(f'Subject Trajectories Across Metrics ({wildcards.c1} vs {wildcards.c2})', fontsize=16)
        
        # Plot each metric with descriptive titles
        metrics_axes = {
            'Tumor_Fraction': 'CNA Tumor Fraction',
            'lt_ctdna_burden': 'FRAGLE Tumor Fraction',
            'mds': 'Motif Diversity Score',
            'nmf': 'Fragment Size NMF',
            'binwise_prob': 'Bin-wise Probability'
        }
        # Add F-profiles to the metrics dictionary
        for fprofile in fprofile_columns:
            metrics_axes[fprofile] = f'F-profile {fprofile.split("_")[1]}'
        
        # Flatten axes array for easier iteration
        axes = axes.flatten()
        
        for i, (metric, title) in enumerate(metrics_axes.items()):
            ax = axes[i]
            # Plot each subject's trajectory
            for subject in subjects:
                subject_data = df[df['subject'] == subject].sort_values('cycle')
                ax.plot(subject_data['cycle'], subject_data[metric], 
                        marker='o', label=subject, color=subject_colors[subject])
            
            ax.set_title(title)
            ax.set_xlabel('Cycle')
            ax.set_ylabel('Value')
            ax.grid(True, alpha=0.3)
            
            # Add legend to the first plot only
            if i == 0:
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Hide any unused subplots
        for i in range(len(all_metrics), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(output.plot, bbox_inches='tight', dpi=300)

rule plot_auc_comparison:
    resources:
        mem_mb = "16000",
        time = "0:10:00",
        cpus = "4"
    input:
        combined_dataset = rules.combine_datasets.output.combined_dataset
    output:
        plot = plot_dir + "/auc_comparison_{c1}_v_{c2}.png"
    log:
        log_dir + "/plot_auc_comparison_{c1}_v_{c2}.log"
    run:
        import matplotlib.pyplot as plt
        import seaborn as sns
        from sklearn.metrics import roc_curve, auc
        # Read the combined dataset
        df = pd.read_csv(input.combined_dataset)
        
        # Filter for the current cohort comparison, excluding healthy controls if they're part of the comparison
        if NORMAL_CONTROL in [wildcards.c1, wildcards.c2]:
            # If healthy is one of the cohorts, only include the other cohort
            other_cohort = wildcards.c2 if wildcards.c1 == NORMAL_CONTROL else wildcards.c1
            df = df[df['cohort'] == other_cohort]
        else:
            # If neither cohort is healthy, include both
            df = df[df['cohort'].isin([wildcards.c1, wildcards.c2])]
        
        # Get all metrics including F-profiles
        base_metrics = ['Tumor_Fraction', 'lt_ctdna_burden', 'mds', 'nmf', 'binwise_prob']
        fprofile_columns = [col for col in df.columns if col.startswith('fprofile_')]
        all_metrics = base_metrics + fprofile_columns
        
        # Convert metrics to numeric, replacing 'missing' with NaN
        for metric in all_metrics:
            df[metric] = pd.to_numeric(df[metric].replace('missing', pd.NA), errors='coerce')
        
        # Calculate AUC for each metric
        auc_scores = {}
        for metric in all_metrics:
            # Create binary labels (1 for c2, 0 for c1)
            labels = (df['cohort'] == wildcards.c2).astype(int)
            # Calculate AUC
            fpr, tpr, _ = roc_curve(labels, df[metric])
            auc_scores[metric] = auc(fpr, tpr)
        
        # Create metric names for display
        metric_names = {
            'Tumor_Fraction': 'CNA Tumor Fraction',
            'lt_ctdna_burden': 'FRAGLE Tumor Fraction',
            'mds': 'Motif Diversity Score',
            'nmf': 'Fragment Size NMF',
            'binwise_prob': 'Bin-wise Probability'
        }
        # Add F-profiles to the names dictionary
        for fprofile in fprofile_columns:
            metric_names[fprofile] = f'F-profile {fprofile.split("_")[1]}'
        
        # Create a DataFrame for plotting
        auc_data = pd.DataFrame({
            'Metric': [metric_names[metric] for metric in all_metrics],
            'AUC': [auc_scores[metric] for metric in all_metrics]
        })
        
        # Sort by AUC in descending order
        auc_data = auc_data.sort_values('AUC', ascending=False)
        
        # Create the bar plot
        plt.figure(figsize=(12, 6))
        ax = sns.barplot(data=auc_data, x='Metric', y='AUC', palette='viridis')
        
        # Add AUC values on top of bars
        for i, row in auc_data.iterrows():
            ax.text(i, row['AUC'] + 0.01, f'{row["AUC"]:.3f}', 
                    ha='center', va='bottom', fontsize=10)
        
        # Customize the plot
        plt.title(f'AUC Comparison Across Metrics ({wildcards.c1} vs {wildcards.c2})', fontsize=14)
        plt.ylim(0, 1.1)  # Add some space above the bars
        plt.ylabel('AUC', fontsize=12)
        plt.xlabel('Metric', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        
        # Add a horizontal line at AUC = 0.5 (random chance)
        plt.axhline(y=0.5, color='r', linestyle='--', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output.plot, dpi=300, bbox_inches='tight')

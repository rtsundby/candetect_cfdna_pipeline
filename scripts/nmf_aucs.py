import argparse
import aucs
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--probs_dir", help='Path to end motifs probs dir (input)', required=True)
parser.add_argument("--max_components", help='Maximum number of components', type=int, required=True)
parser.add_argument("--samplelist", help='Path to samplelist (input)', required=True)
parser.add_argument("--out_dir", help='Path to output directory', required=True)
parser.add_argument("--cohorts", help='Disease cohorts', required=True)

args = parser.parse_args()
print("Running with the following argument values:")
for k,v in sorted(vars(args).items()):
    print(f"\t{k}: {v}")

PROBS_DIR = args.probs_dir
MAX_COMPONENTS = args.max_components
SAMPLELIST = args.samplelist
OUT_DIR = args.out_dir
COHORTS = args.cohorts.split(',')

print(globals().keys())

print(f"Loading dx_df from samplelist: {SAMPLELIST}")
try:
    dx_df = pd.read_csv(SAMPLELIST, sep='\t', header=0)
except Exception as e:
    raise RuntimeError(f"Error loading samplelist file: {e}")
required_columns = ['library', 'cohort'] # Retain only `library` and `cohort` columns
missing_columns = [col for col in required_columns if col not in dx_df.columns]
if missing_columns:
    raise ValueError(f"Missing required columns in samplelist file: {missing_columns}")
dx_df = dx_df[required_columns]
dx_df = dx_df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
print("dx_df:")
print(dx_df.head())

COHORT_1 = []
COHORT_2 = []
for c1, c2 in it.combinations(COHORTS, 2):
    COHORT_1.append(c1)
    COHORT_2.append(c2)
    print(c1, c2)

# create output directories
scores_dir = f"{OUT_DIR}/scores"
aucs_dir = f"{OUT_DIR}/aucs"
aucs_plots_dir = f"{OUT_DIR}/plots/aucs"
aucs_v_comp_plots_dir = f"{OUT_DIR}/plots/aucs_v_components"
os.makedirs(scores_dir, exist_ok=True)
os.makedirs(aucs_dir, exist_ok=True)
os.makedirs(aucs_plots_dir, exist_ok=True)
os.makedirs(aucs_v_comp_plots_dir, exist_ok=True)

cohorts_df = pd.read_csv(SAMPLELIST, sep='\t', index_col=0)

################################################################################
# Compute final score and thresholds for each comparison
# Summarize scores across all comparisons

best_rows = []
rows = []
for c1, c2 in it.combinations(COHORTS, 2):
    best_med_auc = None
    best_row = None
    for n_comp in range(2,MAX_COMPONENTS+1):
        print("c1:", c1, "c2:", c2)
        name = f"{c1}_vs_{c2}.n={n_comp}"
        row = [c1, c2, n_comp]
        
        # get scores df
        probs_df = pd.read_csv(
            f"{PROBS_DIR}/{name}.csv", index_col=0)
        
        # Calculate median score from repeat columns
        repeat_cols = [col for col in probs_df.columns if col.startswith('repeat_')]
        probs_df['med_score'] = probs_df[repeat_cols].median(axis=1)
        
        scores_df = aucs.get_scores_and_thresholds(probs_df,cohorts_df,c1,c2)
        print(scores_df)
        scores_df.to_csv(f"{scores_dir}/{name}.csv")

        # get aucs for each fold, avg_score, med_score
        y_true = scores_df["cohort_int"]
        print(y_true)
        scores_df = scores_df.iloc[:,0:11]
        row += aucs.get_aucs(scores_df, y_true)
        rows.append(row)

        # track best row from this comparison (by med_score auc)
        med_auc = row[-1]
        if (best_med_auc is None) or (med_auc > best_med_auc):
            print("BETTER ROW FOUND:")
            print(row)
            best_med_auc = med_auc
            best_row = row

    best_rows.append(best_row)
    best_med_auc = None
    best_row = None

aucs_df = pd.DataFrame(data=rows,
                    columns=["cohort1", "cohort2","n_components"]+scores_df.columns.to_list())
print(aucs_df)
aucs_df.to_csv(f"{aucs_dir}/aucs.csv")

best_aucs_df = pd.DataFrame(data=best_rows,
                    columns=["cohort1", "cohort2","n_components"]+scores_df.columns.to_list())
print(best_aucs_df)
best_aucs_df.to_csv(f"{aucs_dir}/best_med_aucs.csv")

############################################################################
# Plot AUC summaries
df = best_aucs_df.pivot(index="cohort1", columns="cohort2",
                values="med_score")
df = df.reindex(COHORTS[:-1])
df = df[COHORTS[1:]]
aucs.plot_aucs_heatmap(df,aucs_plots_dir,"med_score")
#df = best_aucs_df.pivot(index="cohort1", columns="cohort2",
#                values="avg_score")
#df = df.reindex(COHORTS[:-1])
#df = df[COHORTS[1:]]
#aucs.plot_aucs_heatmap(df,aucs_plots_dir,"avg_score")

############################################################################
# Plot AUCs per n-component individually
for c1, c2 in it.combinations(COHORTS, 2):
    print("c1:", c1, "c2:", c2)
    name = f"{c1}_vs_{c2}"
    df = aucs_df.loc[(aucs_df["cohort1"] == c1) & (aucs_df["cohort2"] == c2)]
    df = pd.melt(df, id_vars=["n_components"], value_vars=df.columns[-12:])
    df = df.rename(columns={"variable":"repeat","value":"auc"})
    
    repeats_df = df.loc[~df["repeat"].isin(["med_score","avg_score"])]
    med_score_df = df.loc[df["repeat"] == "med_score"]

    # plot scores from repeats
    ax = sns.boxplot(
        data=repeats_df,
        x="n_components",
        y="auc",
        boxprops={'alpha': 0.4},
        fliersize=0,
        dodge=True,
        color='white',
    )
    sns.stripplot(
        data=repeats_df,
        x="n_components",
        y="auc",
        ax=ax,
        jitter=True,
        dodge=True,
        edgecolor='gray',
        color='gray',
        linewidth=1,
        s=4,
        alpha=0.75,
    )
    sns.stripplot(
        data=med_score_df,
        x="n_components",
        y="auc",
        ax=ax,
        edgecolor='gray',
        color='red',
        linewidth=1,
        marker="X",
        s=7,
        alpha=0.75,
    )
    ax.set_ylim(0.0, 1.0)

    plt.title(name)
    plt.tight_layout()
    plt.savefig(f"{aucs_v_comp_plots_dir}/{name}.pdf")
    plt.clf()
    plt.close()

############################################################################
# Plot AUCs per n-component together

# mapping cohorts to corresponding subplot axis coords
c1_mapping = {}
c2_mapping = {}
iter = 0
for c1, c2 in it.combinations(COHORTS, 2):
    c1_mapping[c1]=iter
    c2_mapping[c2]=iter
    iter +=1

fig, axs = plt.subplots(3, 3, figsize=(12, 12))
for c1, c2 in it.combinations(COHORTS, 2):
    ax = axs[c1_mapping[c1], c2_mapping[c2]]
    
    df = aucs_df.loc[(aucs_df["cohort1"] == c1) & (aucs_df["cohort2"] == c2)]
    df = pd.melt(df, id_vars=["n_components"], value_vars=df.columns[-12:])
    df = df.rename(columns={"variable":"repeat","value":"auc"})
    
    repeats_df = df.loc[~df["repeat"].isin(["med_score","avg_score"])]
    med_score_df = df.loc[df["repeat"] == "med_score"]

    # plot scores from repeats
    sns.boxplot(
        data=repeats_df,
        x="n_components",
        y="auc",
        ax=ax,
        boxprops={'alpha': 0.4},
        fliersize=0,
        dodge=True,
        color='white',
    )
    sns.stripplot(
        data=repeats_df,
        x="n_components",
        y="auc",
        ax=ax,
        jitter=True,
        dodge=True,
        edgecolor='gray',
        color='gray',
        linewidth=1,
        s=2,
        alpha=0.75,
    )
    sns.stripplot(
        data=med_score_df,
        x="n_components",
        y="auc",
        ax=ax,
        edgecolor='gray',
        color='red',
        linewidth=1,
        marker="X",
        s=7,
        alpha=0.75,
    )
    ax.set_ylim(0.0, 1.0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xlabel("Number of Components")
    ax.set_ylabel("AUC")

fig.delaxes(axs[1][0])
fig.delaxes(axs[2][0])
fig.delaxes(axs[2][1])
LABEL_SIZE = 30
FONT_WEIGHT = "normal"
fig.text(0.01, 0.83, 'Healthy', va='center',
         rotation='vertical', fontweight=FONT_WEIGHT, fontsize=LABEL_SIZE)
fig.text(0.01, 0.5, 'PN', va='center',
         rotation='vertical', fontweight=FONT_WEIGHT, fontsize=LABEL_SIZE)
fig.text(0.01, 0.2, 'AN', va='center',
         rotation='vertical', fontweight=FONT_WEIGHT, fontsize=LABEL_SIZE)
fig.text(0.17, 0.01, 'PN', ha='center',
         fontweight=FONT_WEIGHT, fontsize=LABEL_SIZE)
fig.text(0.5, 0.01, 'AN', ha='center',
         fontweight=FONT_WEIGHT, fontsize=LABEL_SIZE)
fig.text(0.83, 0.01, 'MPNST', ha='center',
         fontweight=FONT_WEIGHT, fontsize=LABEL_SIZE)

fig.text(0.5, 0.95, 'AUC vs Components', ha='center',
         fontweight=FONT_WEIGHT, fontsize=36)

# fig.suptitle("AUC vs Components", fontsize=24)
plt.tight_layout(pad=5, h_pad=1, w_pad=1)
plt.savefig(f"{aucs_v_comp_plots_dir}/summary.n={MAX_COMPONENTS}.pdf")
plt.clf()
plt.close()

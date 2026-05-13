#!/bin/bash
#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=48g
#SBATCH --time=120:00:00
#SBATCH --no-requeue

module load mamba_install
module load snakemake
source /home/wuy24/bin/myconda
export MAMBA_NO_BANNER=1
#export CONDA_PKGS_DIRS=/data/CCRCCDI/software/env/pkgs

mkdir logs
mkdir stats
#snakemake -s workflows/frag_reads.smk \
#     --configfile config/cpm_nf1.yaml  \
#     --cores 4   \
#     --latency-wait 600 \
#     --use-conda --conda-frontend mamba \
#     -k --stats stats/snakemake_$(date +"%Y%m%d%H%M%S").stats \
#     --rerun-incomplete --restart-times 0 \
#     -j 300 --printshellcmds   \
#     --cluster-config config/cluster.yaml  \
#     --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}"  \
#     >&  logs/snakemake_$(date +"%Y%m%d%H%M%S").log

#snakemake -s workflows/binwise_frags.smk \
#     --configfile config/cpm_nf1.yaml  \
#     --cores 4   \
#     --latency-wait 600 \
#     --use-conda --conda-frontend mamba \
#     -k --stats stats/snakemake_$(date +"%Y%m%d%H%M%S").stats \
#     --rerun-incomplete --restart-times 0 \
#     -j 300 --printshellcmds   \
#     --cluster-config config/cluster.yaml  \
#     --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}"  \
#     >&  logs/snakemake_$(date +"%Y%m%d%H%M%S").log

#snakemake -s workflows/NMF_size.smk \
#     --configfile config/cpm_nf1.yaml  \
#     --cores 4   \
#     --latency-wait 600 \
#     --use-conda --conda-frontend mamba \
#     -k --stats stats/snakemake_$(date +"%Y%m%d%H%M%S").stats \
#     --rerun-incomplete --restart-times 0 \
#     -j 300 --printshellcmds   \
#     --cluster-config config/cluster.yaml  \
#     --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}"  \
#     >&  logs/snakemake_$(date +"%Y%m%d%H%M%S").log

snakemake -s workflows/end_motifs.smk \
     --configfile config/cpm_nf1.yaml  \
     --cores 4   \
     --latency-wait 600 \
     --use-conda --conda-frontend mamba \
     -k --stats stats/snakemake_$(date +"%Y%m%d%H%M%S").stats \
     --rerun-incomplete --restart-times 0 \
     -j 300 --printshellcmds   \
     --cluster-config config/cluster.yaml  \
     --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}"  \
     >&  logs/snakemake_$(date +"%Y%m%d%H%M%S").log


#snakemake -s workflows/ichor.smk \
#     --configfile config/cpm_nf1.yaml  \
#     --cores 4   \
#     --latency-wait 600 \
#     --use-conda --conda-frontend mamba \
#     -k --stats stats/snakemake_$(date +"%Y%m%d%H%M%S").stats \
#     --rerun-incomplete --restart-times 0 \
#     -j 300 --printshellcmds   \
#     --cluster-config config/cluster.yaml  \
#     --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}"  \
#     >&  logs/snakemake_$(date +"%Y%m%d%H%M%S").log

#!/usr/bin/env /usr/local/Anaconda/envs/py3.9/bin/python
import psutil, argparse, os, subprocess, re, sys, itertools, glob
import snakemake
import os
import pandas as pd
from argparse import ArgumentParser
import yaml

#print("Snakemake version:", snakemake.__version__)
# get the path of the script: e.g. /is2/projects/CCR-SF/active/Software/scripts/bin/
bin_path = os.path.dirname(os.path.realpath(__file__))
# get the main script path: e.g. /is2/projects/CCR-SF/active/Software/scripts/
#script_path = os.path.dirname(bin_path)
script_path = bin_path
#print(script_path)
work_path = os.getcwd()
print(work_path)

from argparse import ArgumentParser

#refpath = f"{bin_path}/currentsnake/"
snakepath = f"{bin_path}"
print(f"pipeline source path is: {snakepath}")

class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def main(raw_args=None):
    parser = argparse.ArgumentParser(description="""Help to set up and run the Single Cell Pipelines""", formatter_class=SmartFormatter)

    #parser.add_argument("-wf", "--workflow", metavar="workflow",
    #    action = "store", type=str, required=True, help="snakemake workflow. [Required]")
    parser.add_argument("-cf", "--config", metavar="config",
        action = "store", type=str, required=True, help="config file. [Required]")
    parser.add_argument("-m", "--mode", metavar="mode",
        action = "store", type=str, required=True, 
        help="""pipeline mode for the run
        model for modeling
        predict for prediction""")
    parser.add_argument("-n", "--dryrun", action="store_true", help="dry run")
    parser.add_argument("--rerun", action="store_true",
        help="dry run then prompt for submitting jobs that need to be rerun")
    parser.add_argument("--unlock", action="store_true",
        help="unlock working directory")
    args = parser.parse_args(raw_args)

    if args.mode == 'model':
        workflow = 'workflows/fragmentomics_model.smk'
    elif args.mode == 'predict':
        workflow = 'workflows/fragmentomics_pred.smk'
    elif args.mode == 'align':
        workflow = 'workflows/frag_reads.smk_backup'
    else:
        raise ValueError("Invalid mode. Please choose 'model' or 'predict'.")
    
    configfile = args.config

    submit_script = f'''#!/bin/bash
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

mkdir -p logs
mkdir -p stats

snakemake -s {workflow} \\
    --configfile {configfile} \\
    --cores 4 \\
    --latency-wait 600 \\
    --use-conda --conda-frontend mamba \\
    -k --stats stats/snakemake_$(date +"%Y%m%d%H%M%S").stats \\
    --rerun-incomplete --restart-times 0 \\
    -j 300 --printshellcmds \\
    --cluster-config config/cluster.yaml \\
    --cluster "sbatch --gres {{cluster.gres}} --cpus-per-task {{cluster.threads}} -p {{cluster.partition}} -t {{cluster.time}} --mem {{cluster.mem}} --job-name {{cluster.name}} --output {{cluster.output}} --error {{cluster.error}}" \\
    >& logs/snakemake_$(date +"%Y%m%d%H%M%S").log
'''

    with open("submit.sh", "w") as f:
        f.write(submit_script)

    if args.dryrun:
        snakemake.snakemake(workflow, configfiles=[args.config], unlock=True)
        snakemake.snakemake(workflow, configfiles=[args.config], printshellcmds=True, dryrun=True)  
    elif args.unlock:
        snakemake.snakemake(workflow, configfiles=[args.config], unlock=True)
    elif args.rerun:
        additionalInfo = False
        if additionalInfo:
            print('Exiting now. Fill in additional information in config file, then rerun with ./run_snakemake_sc.py -wf workflow.smk -cf config.yaml --rerun')
        else:
            stdout = sys.stdout
            sys.stdout = open('dryrun_all_commands.txt', 'w')
            snakemake.snakemake(workflow, configfiles=[args.config], printshellcmds=True, dryrun=True, rerun_triggers=["mtime"])
            sys.stdout = stdout
            snakemake.snakemake(workflow, configfiles=[args.config], printshellcmds=True, dryrun=True,rerun_triggers=["mtime"])
            #print('CONFIG FILE:')
            #with open('config/config.py') as f:
            #    print(f.read())

            submit = input('Submit now? (y/n): ')
            while submit.lower() not in ['y', 'n']:
                submit = input('Please answer yes(y) or no(n): ')
            if submit.lower() == 'y':
                subprocess.check_output('sbatch submit.sh', shell=True)
                #print("Submitted snakemake run for %s in %s" % (args.fastq_path, work_path))
                print("Submitted snakemake run in %s" % work_path)

if __name__ == '__main__':
    main()


#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Master workflor for pre-processing and analysis of cfDNA WGS:            #
# Will run pre-processing first (frag_reads.smk) first. If completed w/o error #
#           Will then run end_motifs, binwise, NMF_size in parallel              #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

rule all:
    input:
        "frag_reads.done",
        "end_motifs.done",
        "binwise.done",
        "NMF_size.done"

# Rule to run frag_reads.smk
rule run_frag_reads:
    output:
        "frag_reads.done"
    params:
        configfile="{configfile}"  # Passed from master_workflow.sh
    shell:
        """
        snakemake -s frag_reads.smk --configfile {params.configfile} --profile profile/ && touch {output} || exit 1
        """

# Rules for the three parallel workflows
rule run_end_motifs:
    input:
        "frag_reads.done"
    output:
        "end_motifs.done"
    params:
        configfile="{configfile}"  
    shell:
        """
        snakemake -s end_motifs.smk --nolock --configfile {params.configfile} --profile profile/ && touch {output}
        """

rule run_binwise:
    input:
        "frag_reads.done"
    output:
        "binwise.done"
    params:
        configfile="{configfile}"  
    shell:
        """
        snakemake -s binwise_frags.smk --nolock --configfile {params.configfile} --profile profile/ && touch {output}
        """

rule run_NMF_size:
    input:
        "frag_reads.done"
    output:
        "NMF_size.done"
    params:
        configfile="{configfile}"  
    shell:
        """
        snakemake -s NMF_size.smk --nolock --configfile {params.configfile} --profile profile/ && touch {output}
        """

import os 
import glob

data_dir = config["data_dir"]
size_selected_dir = config["data_dir"] + "/size_selected_bams"
wig_dir = config["data_dir"] + "/wig_files"
ichor_dir = config["data_dir"] + "/analysis/cfdna_wgs/ichorCNA"
log_dir = config["data_dir"] + "/logs/ichor"

os.makedirs(size_selected_dir, exist_ok=True)
os.makedirs(wig_dir, exist_ok=True)
os.makedirs(ichor_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

LIBS = []
for bam in sorted(glob.glob(data_dir + "/bams/*_filt.bam")):
    # Remove the _filt.bam extension to get just the sample name
    sample_name = os.path.basename(bam).replace("_filt.bam", "")
    LIBS.append(sample_name)
print(LIBS)

rule all:
    input:
        expand(size_selected_dir + "/{samples}.ss.bam", samples=LIBS),
        expand(wig_dir + "/{samples}.wig", samples=LIBS),
        expand(ichor_dir + "/{samples}/{samples}.cna.seg", samples=LIBS),
        ichor_dir + "/tf.csv",

rule size_select:
    resources:
        mem_mb = "64000",
        time   = "0:45:00",
        cpus   = "16",
    input:
        bam_file = data_dir + "/bams/{samples}_filt.bam"
    output:
        ss_bam = size_selected_dir + "/{samples}.ss.bam",
    log:
        log_dir + "/size_select/{samples}.log",
    params:
        out_dir = size_selected_dir,
        log_dir = log_dir + "/size_select",
        min_size = config["min_fragment_size"],
        max_size = config["max_fragment_size"],
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}
        
        module load samtools
        samtools view -h {input.bam_file} -@ {resources.cpus} | \
        awk 'substr($0,1,1)=="@" || ($9>= {params.min_size} && $9<={params.max_size}) || ($9<=-{params.min_size} && $9>=-{params.max_size})' | \
        samtools view -b -@ {resources.cpus} > {output.ss_bam}
        samtools index -@ {resources.cpus} {output.ss_bam}
        """

rule read_counter:
    resources:
        mem_mb = "32000",
        time   = "0:25:00",
        cpus   = "16",
    input:
        bam_file = size_selected_dir + "/{samples}.ss.bam",
    output:
        wig_file = wig_dir + "/{samples}.wig",
    log:
        log_dir + "/wig_files/{samples}.log",
    params:
        readCounter = config["readCounter_script"],
        window		= config["window"],
        quality		= config["quality"],
        chrs		= config["chromosomes"],
        out_dir = wig_dir,
        log_dir = log_dir + "/wig_files",
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        {params.readCounter} \
        -w {params.window} \
        -q {params.quality} \
        -c {params.chrs} \
        {input.bam_file} \
        > {output.wig_file} \
        2> {log}
        """
        
rule ichorCNA:
    resources:
        mem_mb = "64000",
        time   = "0:10:00",
        cpus   = "16",
    input:
        tum_wig = wig_dir + "/{tumor}.wig",
    output:
        cna = ichor_dir + "/{tumor}/{tumor}.cna.seg",
        #corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
        #param="results/ichorCNA/{tumor}/{tumor}.params.txt",
        #segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
        #seg="results/ichorCNA/{tumor}/{tumor}.seg",
        #rdata="results/ichorCNA/{tumor}/{tumor}.RData"
    log:
        log_dir + "/ichorCNA/{tumor}.log",
    params:
        rscript           = config["ichorCNA_rscript"],
        libdir            = config["ichorCNA_libdir"], # added
        id                = "{tumor}",
        gcwig             = config["ichorCNA_gcWig"],
        mapwig            = config["ichorCNA_mapWig"],
        centromere        = config["ichorCNA_centromere"],
        normalpanel       = config["ichorCNA_normalPanel"],
        genomeBuild       = config["ichorCNA_genomeBuild"], # added
        genomeStyle       = config["ichorCNA_genomeStyle"], # added
        chrs              = config["ichorCNA_chrs"],
        chrTrain          = config["ichorCNA_chrTrain"],
        estimateNormal    = config["ichorCNA_estimateNormal"],
        estimatePloidy    = config["ichorCNA_estimatePloidy"],
        estimateClonality = config["ichorCNA_estimateClonality"],
        normal            = config["ichorCNA_normal"],
        ploidy            = config["ichorCNA_ploidy"],
        maxCN             = config["ichorCNA_maxCN"],
        scStates          = config["ichorCNA_scStates"],
        txnE              = config["ichorCNA_txnE"],
        txnStrength       = config["ichorCNA_txnStrength"],
        out_dir           = ichor_dir + "/{tumor}",
        log_dir           = log_dir + "/ichorCNA",

        #### UNUSED OPTIONS
        # includeHOMD=config["ichorCNA_includeHOMD"],
        # fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
        # minMapScore=config["ichorCNA_minMapScore"],
        # maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
        # maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
        # exons=config["ichorCNA_exons"],
        # plotFileType=config["ichorCNA_plotFileType"],
        # plotYlim=config["ichorCNA_plotYlim"],
        # 

    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}
        module load R
        
        Rscript {params.rscript} \
        --libdir {params.libdir} \
        --WIG {input.tum_wig} \
        --id {params.id} \
        --gcWig {params.gcwig} \
        --mapWig {params.mapwig} \
        --centromere {params.centromere} \
        --normalPanel {params.normalpanel} \
        --chrs \"{params.chrs}\" \
        --chrTrain \"{params.chrTrain}\" \
        --estimateNormal {params.estimateNormal} \
        --estimatePloidy {params.estimatePloidy} \
        --estimateScPrevalence {params.estimateClonality} \
        --normal \"{params.normal}\" \
        --ploidy \"{params.ploidy}\" \
        --maxCN {params.maxCN} \
        --scStates \"{params.scStates}\" \
        --txnE {params.txnE} \
        --txnStrength {params.txnStrength} \
        --outDir {params.out_dir} \
        --genomeStyle {params.genomeStyle} \
        --genomeBuild {params.genomeBuild} \
        &> {log}
        """

        ### UNUSED OPTIONS:
        # """
        # --includeHOMD {params.includeHOMD} \
        # --exons.bed {params.exons} \
        # --minMapScore {params.minMapScore} \
        # --fracReadsInChrYForMale {params.fracReadsChrYMale} \
        # --maxFracGenomeSubclone {params.maxFracGenomeSubclone} \
        # --maxFracCNASubclone {params.maxFracCNASubclone} \
        # --plotFileType {params.plotFileType} \
        # --plotYLim \"{params.plotYlim}\" \
        # """


rule features:
    resources:
        mem_mb = "64000",
        time   = "0:20:00",
        cpus   = "16",  
    input:
        ichor_dir = directory(ichor_dir),
        cna_files = expand(ichor_dir + "/{samples}/{samples}.cna.seg", samples=LIBS)
    output:
        tf_csv = ichor_dir + "/tf.csv"
    log:
        log_dir + "/features/features.log"
    params:
        script = "scripts/features.py",
        log_dir = log_dir + "/features"
    shell:
        """
        # Create log directory
        mkdir -p {params.log_dir}

        # Check if Python script exists
        if [ ! -f {params.script} ]; then
            echo "ERROR: Python script not found at {params.script}" >> {log} 2>&1
            exit 1
        fi

        # Check input directory exists and has content
        if [ ! -d {input.ichor_dir} ]; then
            echo "ERROR: Input directory {input.ichor_dir} does not exist" >> {log} 2>&1
            exit 1
        fi

        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.tf_csv})

        echo "Starting feature extraction..." >> {log} 2>&1
        echo "Input directory: {input.ichor_dir}" >> {log} 2>&1
        echo "Output file: {output.tf_csv}" >> {log} 2>&1
        
        # Load Python module if needed
        module load python || echo "Warning: Failed to load Python module" >> {log} 2>&1

        # Run the Python script with error handling
        python {params.script} \
            --ichor_dir {input.ichor_dir} \
            --output {output.tf_csv} \
            2>> {log} || {{
            echo "ERROR: Feature extraction failed" >> {log} 2>&1
            exit 1
        }}

        # Verify output was created
        if [ ! -f {output.tf_csv} ]; then
            echo "ERROR: Output file was not created" >> {log} 2>&1
            exit 1
        fi

        echo "Feature extraction completed successfully" >> {log} 2>&1
        """

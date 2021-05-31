#######################################################
# Snakemake pipeline for Multi Sample Variant Calling #
#######################################################

# Minimum version needed for this workflow
from snakemake.utils import min_version
min_version("6.0")

# Import libraries needed for this workflow
import pandas as pd

# Configurations
# --------------
configfile: 'config.yaml'
WORKINGDIR = config['working_dir']
OUTPUTDIR = config['output_dir']

# Binaries
# --------

fastp = config['fastp_binary']
multiqc = config['multiqc_binary']
samtools = config['samtools_binary']
bcftools = config['bcftools_binary']
minimap2 = config['minimap2_binary']

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str, index_col=0, sep="\t")
samplefile = config["units"]

# Functions will be used in the follow through
def is_single_end(sample):
	"""
	This function detect missing value in the column 2 of the units.tsv
	"""
	if "fq2" not in samples.columns:
		return True
	else:
		return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
	"""This function checks if the sample has paired end or single end reads and returns 1 or 2 names of the fastq files"""
	if is_single_end(wildcards.sample):
		return samples.loc[(wildcards.sample), ["fq1"]].dropna()
	else:
		return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trim_names(wildcards):
	"""
	This function:
	  1. Checks if the sample is paired end or single end
	  2. Returns the correct input and output trimmed file names. 
	"""
	if is_single_end(wildcards.sample):
		inFile = samples.loc[(wildcards.sample), ["fq1"]].dropna()
		return "--in1 " + inFile[0] + " --out1 " + WORKINGDIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz" 
	else:
		inFile = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
		return "--in1 " + inFile[0] + " --in2 " + inFile[1] + " --out1 " + WORKINGDIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz --out2 "  + WORKINGDIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"


#################
# Desired outputs
#################
MULTIQC = OUTPUTDIR + "multiqc_report.html"
BAM_FILES = expand(OUTPUTDIR + "bamfiles/{sample}.sort.bam", sample = SAMPLES)
# MAPPING_REPORT = OUTPUTDIR + "mapping_summary.csv"
BCFTOOLS_VCFFILE = OUTPUTDIR + "bcftools.vcf.gz"
DEEPVARIANT_VCFFILE = OUTPUTDIR + "deepvariant.vcf.gz"


rule all:
	input:
		MULTIQC,
		BAM_FILES,
		BCFTOOLS_VCFFILE,
		OUTPUTDIR + "deepvariant.vcf.gz"
	message:
		"Variant calling pipeline run complete..!!"
	shell:
		"cp config.yaml {OUTPUTDIR};"
		"cp samples.tsv {OUTPUTDIR};"


rule reference_symlink:
	input:
		genome_fasta = config['genome_fasta']
	output:
		genome_fasta_symlink = WORKINGDIR + 'reference/genome.fasta'
	params:
		genome_dir = WORKINGDIR + 'reference/'
	message:
		"Creating reference genome symlink..!!"
	shell:
		# "mkdir {params.genome_dir}; "
		"ln -s {input.genome_fasta} {output.genome_fasta_symlink}"

rule minimap_index:
	input:
		genome_fasta = WORKINGDIR + 'reference/genome.fasta'
	output:
		genome_index = WORKINGDIR + 'reference/genome.mmi'
	message:
		"Creating minimap2 index for the reference genome..!!"
	log:
		"logs/minimap2/index/stdout.log"
	shell:
		"{minimap2} -d {output.genome_index} {input.genome_fasta} 2> {log}"

rule samtools_index:
	input:
		WORKINGDIR + 'reference/genome.fasta'
	output:
		WORKINGDIR + 'reference/genome.fasta.fai'
	message:
		"Creating samtools fai index for the reference genome..!!"
	log:
		"logs/samtools/reference_index.log"
	wrapper:
		"0.74.0/bio/samtools/faidx"


rule fastp:
	input:
		get_fastq
	output:
		fq1  = temp(WORKINGDIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz"),
		fq2  = temp(WORKINGDIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"),
		html = WORKINGDIR + "fastp/{sample}_fastp.html",
		json = WORKINGDIR + "fastp/{sample}_fastp.json"
	message:
		"trimming {wildcards.sample} reads"
	threads:
		10
	message:
		"fastp trimming for the {sample}..!!"
	log:
		"logs/fastp/{sample}.log"
	params:
		sampleName = "{sample}",
		in_and_out_files =  get_trim_names,
		qualified_quality_phred = config['fastp']['qualified_quality_phred'],
		adaptor_fasta = config['fastp']['adaptor_fasta']
	resources:
		cpus = 10
	run:
		if is_single_end(wildcards.sample) :
			shell(
				"touch {output.fq2}; {fastp} --thread {threads} \
				--qualified_quality_phred 20 \
				--unqualified_percent_limit 30 \
				--cut_tail_window_size 1 \
				--cut_tail_mean_quality 3 \
				--cut_tail \
				--cut_front \
				--cut_front_mean_quality 3 \
				--cut_front_window_size 1 \
				--length_required 35 \
				--html {output.html} \
				--json {output.json} \
				{params.in_and_out_files} \
				2> {log}")
		else:
			shell(
				"{fastp} --thread {threads} \
				--adapter_fasta {params.adaptor_fasta} \
				--correction \
				--qualified_quality_phred 20 \
				--unqualified_percent_limit 30 \
				--cut_tail_window_size 1 \
				--cut_tail_mean_quality 3 \
				--cut_tail \
				--cut_front \
				--cut_front_mean_quality 3 \
				--cut_front_window_size 1 \
				--length_required 35 \
				--html {output.html} \
				--json {output.json} \
				{params.in_and_out_files} \
				2> {log}")

rule multiqc:
	input:
		expand(WORKINGDIR + "fastp/{sample}_fastp.json", sample = SAMPLES)
	output:
		OUTPUTDIR + "multiqc_report.html"
	params:
		# fastp_directory = WORKINGDIR + "fastp/",
		# fastp_directory = lambda w, input: os.path.splitext(input[0])[0],
		# outdir = OUTPUTDIR
	message:
		"Summarising fastp reports with multiqc..!!"
	# shell:
	# 	"{multiqc} --force "
	# 	"--outdir {params.outdir} "
	# 	"{params.fastp_directory} "
	log:
		"logs/multiqc.log"
	wrapper:
		"0.74.0/bio/multiqc"

rule minimap_mapping:
	input:
		genome_index = WORKINGDIR + "reference/genome.mmi",
		genome_fasta = WORKINGDIR + "reference/genome.fasta",
		read1 = WORKINGDIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
		read2 = WORKINGDIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"
	output:
		temp(OUTPUTDIR + "bamfiles/" + "{sample}.sort.bam")
	params:
		readgroup = "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'"
	threads:
		10
	message:
		"Miminap2 mapping of {wildcards.sample} on the reference genome..!!"
	log:
		"logs/minimap2/mapping/{sample}.log"
	run:
		if is_single_end(wildcards.sample):
			shell(
			"{minimap2} -ax sr \
			-R {params.readgroup} \
			-t {threads} \
			-2 \
			--secondary=no \
			{input.genome_index} \
			{input.read1} 2> {log}| {samtools} sort \
			-@ {threads} \
			--reference {input.genome_fasta} \
			-o {output} - 2> {log}")
		else:
			shell(
			"{minimap2} -ax sr \
			-R {params.readgroup} \
			-t {threads} \
			-2 \
			--secondary=no \
			{input.genome_index} \
			{input.read1} \
			{input.read2} 2> {log} | {samtools} sort \
			-@ {threads} \
			--reference {input.genome_fasta} \
			-o {output} - 2> {log}")

rule list_bams:
	input:
		expand(OUTPUTDIR + "bamfiles/" + "{sample}.dedup.bam", sample=SAMPLES),
		expand(OUTPUTDIR + "bamfiles/" + "{sample}.dedup.bam.bai", sample=SAMPLES)
	output:
		temp(OUTPUTDIR + "SortedBams")
	params:
		bam_dir = OUTPUTDIR + "bamfiles/*.dedup.bam"
	message:
		"Creating bam list file..!!"
	shell:
		r"ls -c1 {params.bam_dir} > {output}"

rule bam_index:
	input:
		OUTPUTDIR + "bamfiles/{sample}.sort.bam"
	output:
		temp(OUTPUTDIR + "bamfiles/{sample}.sort.bam.bai")
	threads:
		8
	message:
		"Building bam index for the sample {wildcards.sample}..!!"
	log:
		"logs/samtools/bamindex/{sample}.log"
	wrapper:
		"0.74.0/bio/samtools/index"

rule mark_duplicates:
    input:
        OUTPUTDIR + "bamfiles/{sample}.sort.bam"
    output:
        bam = OUTPUTDIR + "bamfiles/{sample}.dedup.bam",
        metrics = "dedup/{sample}.metrics.txt"
    params:
        "REMOVE_DUPLICATES=true"
    resources:
        mem_mb=100000
    message:
    	"Marking duplicates for the sample {wildcards.sample}..!!"
    log:
        "logs/picard/dedup/{sample}.log"
    wrapper:
        "0.74.0/bio/picard/markduplicates"

rule dedup_bam_index:
	input:
		OUTPUTDIR + "bamfiles/{sample}.dedup.bam"
	output:
		OUTPUTDIR + "bamfiles/{sample}.dedup.bam.bai"
	threads:
		8
	message:
		"Building bam index for the sample {wildcards.sample}..!!"
	log:
		"logs/samtools/bamindex/{sample}.log"
	wrapper:
		"0.74.0/bio/samtools/index"

rule call_vcf:
	input:
		genome_fasta = WORKINGDIR + "reference/genome.fasta",
		sortbamlist = OUTPUTDIR + "SortedBams"
	output:
		BCFTOOLS_VCFFILE
	params:
		"" # optional params string
	message:
		"Calling variants using bcftools..!!"
	log:
		"logs/bcftools/mpileup/stdout.log"
	threads:
		workflow.cores * 0.5
	shell:
		"{bcftools} mpileup -a DP,AD -O u --threads {threads} -f {input.genome_fasta} --bam-list {input.sortbamlist} 2> {log} | {bcftools} call --threads {threads} -mv -O z -o {output} 2>> {log}"

rule deepvariant_gvcf:
	input:
		bam = OUTPUTDIR + "bamfiles/{sample}.dedup.bam",
		bam_index = OUTPUTDIR + "bamfiles/{sample}.dedup.bam.bai",
		ref = WORKINGDIR + "reference/genome.fasta",
		ref_index = WORKINGDIR + "reference/genome.fasta.fai"
	output:
		vcf = OUTPUTDIR + "gvcf_calls/{sample}.vcf.gz",
		gvcf = OUTPUTDIR + "gvcf_calls/{sample}.g.vcf.gz"
	params:
		model="wgs",   # {wgs, wes, pacbio, hybrid}
		sample_name = lambda w: w.sample,
		extra=""
	threads: 
		workflow.cores * 0.5
	message:
		"Calling gvcf for the sample {wildcards.sample} using deepvariant..!!"
	log:
		"logs/deepvariant/{sample}/stdout.log"
	wrapper:
		"0.74.0/bio/deepvariant"

rule glnexus:
	input:
		gvcfs=lambda w: expand(OUTPUTDIR + "gvcf_calls/{sample}.g.vcf.gz", sample=SAMPLES),
	output:
		vcf = DEEPVARIANT_VCFFILE,
		scratch=temp(
			directory("results/all_group_samples_joint_calls/deepvariant.DB")
		),
	threads:
		10
	message:
		"Merging gvcf to vcf using GLnexus..!!"
	container:
		"docker://quay.io/mlin/glnexus:v1.3.1"
	log:
		"logs/glnexus/stdout.log"
	shell:
		"glnexus_cli "
		"--config DeepVariant "
		"--dir {output.scratch} "
		"--threads {threads} "
		"{input.gvcfs} "
		"2> {log} "
		"| bcftools view - "
		"| bgzip -c "
		"> {output.vcf}"

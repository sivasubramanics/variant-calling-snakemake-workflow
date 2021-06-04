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
configfile: 'config/config.yaml'
WORKINGDIR = config['working_dir']
OUTPUTDIR = config['output_dir']

# Binaries
# --------
fastp = config['fastp_binary']
samtools = config['samtools_binary']
bcftools = config['bcftools_binary']
minimap2 = config['minimap2_binary']

# read the tab separated table containing the sample and fastq file names	
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names 
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str, index_col=0, sep="\t")
samplefile = config["units"]


###################
# Desired outputs #
###################
MULTIQC = OUTPUTDIR + "multiqc_report.html"
BAM_FILES = expand(OUTPUTDIR + "bamfiles/{sample}.sort.bam", sample = SAMPLES)
BCFTOOLS_VCFFILE = OUTPUTDIR + "bcftools.vcf.gz"
DEEPVARIANT_VCFFILE = OUTPUTDIR + "deepvariant.vcf.gz"



# Rule defining expected outputfiles.
rule all:
	input:
		MULTIQC,
		BAM_FILES,
		BCFTOOLS_VCFFILE,
		DEEPVARIANT_VCFFILE
	message:
		"Variant calling pipeline run complete..!!"
	shell:
		"cp config.yaml {OUTPUTDIR};"
		"cp samples.tsv {OUTPUTDIR};"

include: "rules/common.smk"
include: "rules/ref.smk"
include: "rules/readprocessing.smk"
include: "rules/mapping.smk"
include: "rules/varcall-bcftools.smk"
include: "rules/varcall-deepvariant.smk"



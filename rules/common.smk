import pandas as pd


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



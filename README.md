# Snakemake Workflow: Variant Calling using bcftools and DeepVaraint

This is a Snakemake workflow implementing the multi-sample variant calling with "bcftools" and DeepVariant. It also allows for single-sample variant calling, and creating a merged vcf files.

## Usage

This workflow requires few preinstalled tools, which help few rules handle pipeing of STDOUT in order to reduce processing time.

### Install these prerequisites:

#### Option 1 (Recommended)

Creating a conda environment to handle the dependecies smooth using [mamba](https://mamba.readthedocs.io/en/latest/installation.html)

```bash
 mamba create -n snakemake -c conda-forge -c bioconda snakemake python=3.6 singularity fastp minimap2 samtools bcftools
```
Note: There maybe dependencies failure if we install these tools using pip or apt due to the version mismatch

#### Option 2

- [x] snakemake: [min version 6.0](https://snakemake.readthedocs.io/en/v6.0.0/getting_started/installation.html#:~:text=for%20installing%20Snakemake.-,Installation%20via%20Conda,See%20here%20for%20installation%20instructions.)
- [x] fastp: [Installation reference url](https://github.com/OpenGene/fastp)
- [x] minimap2: [Installation reference url](https://github.com/lh3/minimap2)
- [x] samtools: [Installation reference url](https://github.com/samtools/samtools)
- [x] bcftools: [Installation reference url](https://github.com/samtools/bcftools)


### config.yaml
Update the path of above mentioned executables in the config/config.yaml as below:

```bash
fastp_binary: /data/Programs/fastp/fastp
minimap2_binary: /data/Programs/minimap2/minimap2
samtools_binary: /data/Programs/samtools-1.9/samtools
bcftools_binary: /data/Programs/bcftools-1.9/bcftools
```

Modifye the config.yaml for the referene genome fasta file.

```bash
genome_fasta: test-data/NC_045512.fasta
```

Further you can change the output and/or working directory path provided in the config.yaml file.

### samples.tsv
Update the config/samples.tsv

```bash
sample	fq1
sample1	test-data/sample1.fq.gz
sample2	test-data/sample2.fq.gz
```
sample| fq1
------|----
sample1|test-data/sample1.fq.gz
sample2|test-data/sample2.fq.gz

Try running the Snakemake workflow for dry run:
```bash
snakemake -npr
```
If there is no error found in the above step, probably we are good to go.

Finally the actual run:
```bash
snakemake --cores all --use-conda --use-singularity
```
This above line will initiate the workflow and generate the outputs and log files.

### Contact:

[@cssivasubramani](https://github.com/sivasubramanics) | c.s.sivasubramani@gmail.com


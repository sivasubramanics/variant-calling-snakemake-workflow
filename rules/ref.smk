import os


# Rule to make symbolic link for the reference genome
rule reference_symlink:
	input:
		genome_fasta = os.path.abspath(config['genome_fasta'])
	output:
		genome_fasta_symlink = WORKINGDIR + 'reference/genome.fasta'
	params:
		genome_dir = WORKINGDIR + 'reference/'
	message:
		"Creating reference genome symlink..!!"
	shell:
		# "mkdir {params.genome_dir}; "
		"ln -s {input.genome_fasta} {output.genome_fasta_symlink}"

# Rule to index reference genome using minimap2
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


# Rule to index reference genome using samtools.
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

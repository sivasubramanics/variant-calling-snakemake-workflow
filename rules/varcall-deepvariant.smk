
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

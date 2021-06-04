# Rule to do QC and trim of fastq reads using fastp.
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

# Rule to compile fastp results.
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

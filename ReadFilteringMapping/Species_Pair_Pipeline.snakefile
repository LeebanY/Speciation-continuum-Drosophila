SAMPLES, = glob_wildcards("{sample}_1.fastq")
configfile: "snakemake_dros.yml"
rule all:
	input:
		expand("addreadgrp_reads/{sample}.bam.bai", sample=SAMPLES)

ruleorder: fastp_pe > single_end_fastp
ruleorder: bwa_map_pe > bwa_map_se

rule single_end_fastp:
	input:
		in_read1="{sample}_1.fastq"
	output:
		out_read1="cleaned/{sample}_1.cleaned.fastq"
	threads: 16
	log:
		"single_end_pe.{sample}.log"
	shell:
		"fastp -i {input.in_read1} -o {output.out_read1}"

rule fastp_pe:
	input:
		in_read1="{sample}_1.fastq", in_read2="{sample}_2.fastq"
	output:
		out_read1="cleaned/{sample}_1.cleaned.fastq", out_read2="cleaned/{sample}_2.cleaned.fastq"
	threads: 16
	log:
		"fast_pe.{sample}.log"
	shell:
		"fastp -i {input.in_read1} -I {input.in_read2} -o {output.out_read1} -O {output.out_read2}"

rule bwa_map_se:
	input:
		"cleaned/{sample}_1.cleaned.fastq"
	output:
		"mapped_reads/{sample}.sam"
	threads: 16
	log:
		"bwa_map_se.{sample}.log"
	shell:
		"bwa mem -x ont2d -t 16 refgenB {input} > {output}"

rule bwa_map_pe:
	input:
		# Need to change reference genome for each independent run
		"cleaned/{sample}_1.cleaned.fastq","cleaned/{sample}_2.cleaned.fastq" 
	output:
		"mapped_reads/{sample}.sam"
	threads: 16
	log:
		"bwa_map_pe.{sample}.log"
	shell:
		"bwa-mem2 mem -p refgenA -t 16 {input} > {output}"

rule convert_bam:
	input:
		"mapped_reads/{sample}.sam"
	output:
		"mapped_reads/{sample}.bam"
	threads: 16
	shell:
		"sambamba view -t 16 -S -f bam {input} > {output}"


rule sambamba_sort:
	input:
		"mapped_reads/{sample}.bam"
	output:
		"sorted_reads/{sample}.bam"
	threads:16
	shell:
		"sambamba sort -t 16 -o {output} {input}"



rule sambamba_markdup:
	input:
		"sorted_reads/{sample}.bam"
	output:
		"markdup_reads/{sample}.bam"
	threads:16
	shell:
		"sambamba markdup -t 16 -r {input} {output}"


rule flagstat:
	input:
		"markdup_reads/{sample}.bam"
	output:
		"flagstat/{sample}.txt"
	shell:
		"sambamba flagstat -t 16 {input} > {output}"

rule sambamba_index:
	input:
		"markdup_reads/{sample}.bam"
	output:
		"markdup_reads/{sample}.bam.bai"
	threads:16
	shell:
		"sambamba index -t 16 {input} {output}"

rule addreadgrps:
	input:
		"markdup_reads/{sample}.bam"
	output:
		"addreadgrp_reads/{sample}.bam"
	params:
		"RGID={sample} RGLB=dros RGPL=illumina RGPU=unit1 RGSM={sample}"
	log:
		"{sample}.readgrp.log"
	wrapper:
		"v0.80.1/bio/picard/addorreplacereadgroups"

rule readgrp_index:
	input:
		"addreadgrp_reads/{sample}.bam"
	output:
		"addreadgrp_reads/{sample}.bam.bai"
	threads:16
	shell:
		"sambamba index -t 16 {input} {output}"

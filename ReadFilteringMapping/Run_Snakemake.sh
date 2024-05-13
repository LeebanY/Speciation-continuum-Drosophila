#!/bin/bash
#SBATCH --job-name=speccontin    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=10G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=specconti.log   # Standard output and error log
pwd; hostname; date

#ACTIVATE CORRECT ENVIRONMENT.
conda activate snakemake_drosophilomics

#REFERENCE should the only one with .fas as ending.
ref=$(ls *fas)

suffix=".fas";
ref_suf=${ref%$suffix}

bwa-mem2 index -p $ref_suf $ref
#bwa index $ref

# MAKE SURE CORRECT REFERENCE IS USED.
sed -i "s/refgenA/$ref_suf/g" Species_Pair_Pipeline.snakefile
#sed -i "s/refgenB/$ref/g" Species_Pair_Pipeline.snakefile

mkdir addreadgrp_reads

# CALL SNAKEMAKE -- OUTPUT BAM FILE WITH READ GROUPS AS SAMPLE NAME.
snakemake -s Species_Pair_Pipeline.snakefile --cluster "sbatch --job-name=snakemake_test --nodes=1 --ntasks=16 --partition=medium --mem=20G" -j 1 --conda-frontend conda --until addreadgrps

mkdir calls

# OUTPUT BAM LIST.
ls *.fastq | cut -d'_' -f1 | uniq | sed 's|^|addreadgrp_reads/|g' | sed 's/$/.bam/g' > bamlist.txt

# VCF OUTPUT NAME.
vcf_output=$(ls *.fastq | cut -d'_' -f1 | uniq | tr '\n' '_' | sed 's/.$//')

#INDEX REF
samtools faidx $ref

# RUN FREEBAYES.
freebayes -f $ref --haplotype-length -1 --no-population-priors \
--hwe-priors-off --use-mapping-quality --ploidy 2 --theta 0.02 --bam-list bamlist.txt > calls/calls.combined.vcf

# ZIP AND INDEX VCF.
bgzip -c calls/calls.combined.vcf > calls/calls.combined.vcf.gz
tabix -p vcf calls/calls.combined.vcf.gz

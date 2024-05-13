#!/bin/bash
#SBATCH --job-name=speccontin    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=10G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=specconti.log   # Standard output and error log
pwd; hostname; date

#conda activate snakemake_drosophilomics
#REFERENCE should the only one with .fas as ending.
#ref=$(ls *fas)

#suffix=".fas";
#ref_suf=${ref%$suffix}

#bwa-mem2 index -p $ref_suf $ref
#bwa index $ref

#sed -i "s/refgenA/$ref_suf/g" Species_Pair_Pipeline.snakefile
#sed -i "s/refgenB/$ref/g" Species_Pair_Pipeline.snakefile

#mkdir addreadgrp_reads

#snakemake -s Species_Pair_Pipeline.snakefile --cluster "sbatch --job-name=snakemake_test --nodes=1 --ntasks=16 --partition=medium --mem=20G" -j 1 --conda-frontend conda --until addreadgrps

#mkdir calls

#ls *.fastq | cut -d'_' -f1 | uniq | sed 's|^|addreadgrp_reads/|g' | sed 's/$/.bam/g' > bamlist.txt

#vcf_output=$(ls *.fastq | cut -d'_' -f1 | uniq | tr '\n' '_' | sed 's/.$//')

#samtools faidx $ref

#freebayes -f $ref --haplotype-length -1 --no-population-priors \
#--hwe-priors-off --use-mapping-quality --ploidy 2 --theta 0.02 --bam-list bamlist.txt > calls/calls.combined.vcf

#bgzip -c calls/calls.combined.vcf > calls/calls.combined.vcf.gz
#tabix -p vcf calls/calls.combined.vcf.gz



############################# -- GIMBLE script -----###########################



conda activate gimble

cd calls/


ref=../*.fas
BAM_DIR=../addreadgrp_reads/
VCF=calls.combined.vcf
#Preprocess vcf, bam and reference.
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/PILOT_PAIRS/GIMBLE/gimble/./gIMble preprocess -f $ref -v $VCF -b $BAM_DIR -t 4 -o Gimbled_pair_MINDEPTH4

#Set up the gimbled file.
#Before doing so, make the samples file.
#ls ../*.fastq | cut -d'_' -f1 | uniq | sed 's/$/,/' | awk '{print $1,$1}'| sed 's/.$//' | sed 's/ //g' | sed 's|../||g' > Gimbled_pair.manual.samples.csv
#Now setup the gimble folder for all the information to be stored
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/PILOT_PAIRS/GIMBLE/gimble/./gIMble setup -v Gimbled_pair.vcf.gz -b Gimbled_pair.bed -g Gimbled_pair.genomefile -s Gimbled_pair.manual.samples.csv -o Gimbled_pair_analysis_MINDEPTH4 -f
#Cut blocks
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/PILOT_PAIRS/GIMBLE/gimble/./gIMble blocks -z Gimbled_pair_analysis_MINDEPTH4.z/ -l 200 -M 20000
#Produce summary stats information for the whole genome.
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/PILOT_PAIRS/GIMBLE/gimble/./gIMble info -z Gimbled_pair_analysis_MINDEPTH4..z/ > Summarystats.wholegenome.MINDEPTH4..info.txt


######################### ---- Set up partition for intron --------###########

intron_bed=../*.bed

#Filter intersected bed by only keeping positions that are found in both samples.
awk '$4==2' Gimbled_pair.bed > Gimbled_pair.filtered.bed
#next slop 10bp off each enty in the intron bed 
bedtools slop -i $intron_bed -g Gimbled_pair.genomefile -b -10 > Slopped.bed
#intersect intron bed with gimbled bed 
bedtools intersect -a Slopped.bed -b Gimbled_pair.filtered.bed > Intronic_regions.callable.slopped.bed
#Next we run bed_pruner
python /home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/ACTUAL_PAIRS/GENOME_BED_INFO/bed_pruner.v1.py -b Intronic_regions.callable.slopped.bed -s 0 -m 200 -M 20000 -l 200
#Output of this a file called: Intronic_regions.callable.slopped.l_200.m_200.M_20000.s_0.r_12345.pruned.bed
sort -Vk1 -k2,3 Intronic_regions.callable.slopped.l_200.m_200.M_20000.s_0.r_12345.pruned.bed | bedtools merge -i - > Intronic_regions.callable.slopped.l_200.m_200.M_20000.s_0.r_12345.pruned.sorted.merged.bed
#Then we intersect callable regions with pruned intronic regions that are callable.
bedtools intersect -a Gimbled_pair.filtered.bed -b Intronic_regions.callable.slopped.l_200.m_200.M_20000.s_0.r_12345.pruned.sorted.merged.bed > Gimbled.callable.intronic.segments.bed
#You need to filter the genomefile to only include those contigs that actually remain.
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' Gimbled.callable.intronic.segments.bed Gimbled_pair.genomefile > Gimbled_pair.genomefile.filtered
#Now set up the gimble analysis again
#Bedfile to use: Gimbled.callable.intronic.segments.bed
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/PILOT_PAIRS/GIMBLE/gimble/./gIMble setup -v Gimbled_pair.vcf.gz -b Gimbled.callable.intronic.segments.bed -g Gimbled_pair.genomefile.filtered -s Gimbled_pair.manual.samples.csv -o Gimbled_pair_analysis -f
#Cut blocks
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/PILOT_PAIRS/GIMBLE/gimble/./gIMble blocks -z Gimbled_pair_analysis.z/ -l 200 -m 20000
#Produce summary stats information for the whole genome.
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/PILOT_PAIRS/GIMBLE/gimble/./gIMble info -z Gimbled_pair_analysis.z/ > Summarystats.introns.info.txt
#Once blocks have been cut export the blocks to tsv file for demographic analysis
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/PILOT_PAIRS/GIMBLE/gimble/./gIMble query -z Gimbled_pair_analysis.z/ -b --bsfs

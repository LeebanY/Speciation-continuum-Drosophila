#!/bin/bash
#SBATCH --job-name=exon_check    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=10G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=dataforallpartitions_deptfilt2.log   # Standard output and error log
pwd; hostname; date

conda activate gimble

cd calls/
ref=../*.fas
samples_file=Gimbled_pair.manual.samples.csv
#intron_bed_FILT=Gimbled.FINALMIND2.callable.intronic.segments.bed
#intron_gfile=Gimbled_pair_FINALMIND2.genomefile.filtered

#PREPROCESS STEP -- REDO. This was the critical step that resulted in ireggular up outputs
#/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble preprocess -f ../*.fas -v calls.combined.vcf.gz -b ../addreadgrp_reads/ -m 2 -o FINAL_dir_mindp2 -t 16
# WHOLE GENOME
#/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble setup -v FINAL_dir_mindp2.vcf.gz -b FINAL_dir_mindp2.bed -g FINAL_dir_mindp2.genomefile -s $samples_file -o $TMPDIR/FINAL_dir_mindp2_WG -f
#Cut blocks
#/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble blocks -z $TMPDIR/FINAL_dir_mindp2_WG.z/ -l 200 -m 20000
#Produce summary stats information for the whole genome.
#/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble info -z $TMPDIR/FINAL_dir_mindp2_WG.z/ > Summarystats.wholegenome.FINALMINDP2.info.txt

echo "-----------Completed processing of vcf and whole genome statistics.-----------------"

#INTRONS -- MIN DP2

#Now set up the gimble analysis again
intron_bed=../*introns.bed

#Filter intersected bed by only keeping positions that are found in both samples.
#awk '$4==2' FINAL_dir_mindp2.bed > FINAL_dir_mindp2.filtered.bed
#next slop 10bp off each enty in the intron bed
#bedtools slop -i $intron_bed -g FINAL_dir_mindp2.genomefile -b -10 > Slopped_FINALMINDP2.bed
#intersect intron bed with gimbled bed
#bedtools intersect -a Slopped_FINALMINDP2.bed -b FINAL_dir_mindp2.filtered.bed > Intronic_regions.FINALMIND2.callable.slopped.bed
#Next we run bed_pruner
python /home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/ACTUAL_PAIRS/GENOME_BED_INFO/bed_pruner.v1.py -b Intronic_regions.FINALMIND2.callable.slopped.bed -s 0 -m PLACEHOLDER -M 20000 -l PLACEHOLDER
#Sort the output of bed_pruner.
sort -Vk1 -k2,3 Intronic_regions.FINALMIND2.callable.slopped.l_PLACEHOLDER.m_PLACEHOLDER.M_20000.s_0.r_12345.pruned.bed | bedtools merge -i - > Intronic_regions.FINALMIND2.callable.slopped.l_PLACEHOLDER.m_PLACEHOLDER.M_20000.s_0.r_12345.pruned.sorted.merged.bed
#Then we intersect callable regions with pruned intronic regions that are callable.
bedtools intersect -a FINAL_dir_mindp2.filtered.bed -b Intronic_regions.FINALMIND2.callable.slopped.l_PLACEHOLDER.m_PLACEHOLDER.M_20000.s_0.r_12345.pruned.sorted.merged.bed > Gimbled.FINALMIND2.callable.intronic.segments.bed
#You need to filter the genomefile to only include those contigs that actually remain.
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' Gimbled.FINALMIND2.callable.intronic.segments.bed FINAL_dir_mindp2.genomefile > Gimbled_pair_FINALMIND2.genomefile.filtered

echo "Completed processing of intermediate intron files and genomefile."

bedtools intersect -a Gimbled.FINALMIND2.callable.intronic.segments.bed -b INTRONS_1_INTERVAL_PER_GENE.bed > Gimbled.FINALMIND2.callable.intronic.segments.FINAL.bed
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' Gimbled.FINALMIND2.callable.intronic.segments.FINAL.bed FINAL_dir_mindp2.genomefile > Gimbled_pair_FINALMIND2.genomefile.FINAL.filtered

intron_bed_FILT=Gimbled.FINALMIND2.callable.intronic.segments.FINAL.bed
intron_gfile=Gimbled_pair_FINALMIND2.genomefile.FINAL.filtered

#Now set up the gimble analysis again
#Bedfile to use: Gimbled.callable.intronic.segments.bed
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble setup -v FINAL_dir_mindp2.vcf.gz -b $intron_bed_FILT -g $intron_gfile -s Gimbled_pair.manual.samples.csv -o $TMPDIR/FINAL_dir_mindp2_INTRONS -f
#Cut blocks
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble blocks -z $TMPDIR/FINAL_dir_mindp2_INTRONS.z/ -l PLACEHOLDER -m 20000
#Produce summary stats information for the whole genome.
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble info -z $TMPDIR/FINAL_dir_mindp2_INTRONS.z/ > Summarystats.introns.FINALMINDP2.furtherfilt.info.txt
#Produce the bsfs
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble query -z $TMPDIR/FINAL_dir_mindp2_INTRONS.z/ -b --bsfs

echo "END SCRIPT"

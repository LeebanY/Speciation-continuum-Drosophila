#!/bin/bash
#SBATCH --job-name=exon_check    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=10G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=dataforallpartitions_deptfilt2.log   # Standard output and error log
pwd; hostname; date

conda activate gimble

ref=../*.fas
samples_file=Gimbled_pair.manual.samples.csv

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

echo "END SCRIPT"

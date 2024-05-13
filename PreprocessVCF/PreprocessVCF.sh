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

#PREPROCESS STEP -- REDO. This was the critical step that resulted in ireggular up outputs
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble preprocess -f ../*.fas -v calls.combined.vcf.gz -b ../addreadgrp_reads/ -m 2 -o FINAL_dir_mindp2 -t 16
# WHOLE GENOME
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble setup -v FINAL_dir_mindp2.vcf.gz -b FINAL_dir_mindp2.bed -g FINAL_dir_mindp2.genomefile -s $samples_file -o $TMPDIR/FINAL_dir_mindp2_WG -f
#Cut blocks
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble blocks -z $TMPDIR/FINAL_dir_mindp2_WG.z/ -l 200 -m 20000
#Produce summary stats information for the whole genome.
/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/GIMBLE/gIMble/./gIMble info -z $TMPDIR/FINAL_dir_mindp2_WG.z/ > Summarystats.wholegenome.FINALMINDP2.info.txt

echo "-----------Completed processing of vcf and whole genome statistics.-----------------"
echo "END SCRIPT"

#!/bin/bash
#SBATCH --job-name=outputdist    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=10G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=outputdist.log   # Standard output and error log
pwd; hostname; date

conda activate gimble

cd calls/
ref=../*.fas
samples_file=Gimbled_pair.manual.samples.csv

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

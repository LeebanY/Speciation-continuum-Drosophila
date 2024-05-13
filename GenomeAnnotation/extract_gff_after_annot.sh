#!/bin/bash
#SBATCH --job-name=tidy_gff    # Job name
#SBATCH --ntasks=4                    # Run on a single CPU
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=tidygff.log   # Standard output and error log
pwd; hostname; date

conda activate genometools

#D.bipectinata.fasta.fixedheaders.fasta.masked_braker 
#File to loop through all braker directories and retrieve the gff file. The braker file should then be reworked a bit to include introns.
#These files should then be processed to get out slopped introns ready for the analysis. 

for all in BRAKER_annotations_species/*;do
	fbname=$(basename "$all" | cut -d'.' -f1-3)
	#echo $fbname
	#gt gff3 -sort -fixregionboundaries -checkids -tidy -setsource BRAKER $all/augustus.hints.gff3 | perl -pe 's/_\[organism\S+//g' > Species_Annotations/$fbname.braker.gt.gff3
	#gff2bed < Species_Annotations/$fbname.braker.gt.gff3 > Species_Annotations/$fbname.braker.gt.bed
	awk '$8=="intron"' Species_Annotations/$fbname.braker.gt.bed > Species_Annotations/introns/$fbname.braker.gt.introns.bed
done

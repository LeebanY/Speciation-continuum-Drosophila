#!/bin/bash
#SBATCH --job-name=BUSCO    # Job name
#SBATCH --array=1-8
#SBATCH --ntasks=6                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=BUSCO.log   # Standard output and error log
pwd; hostname; date

conda activate busco

FILES=(renamed_headers_fasta/*.fas)
busco -i ${FILES[$SLURM_ARRAY_TASK_ID]} -l diptera_odb10 -o ${FILES[$SLURM_ARRAY_TASK_ID]}.BUSCO -m genome -c 6 -f

#busco -i D.ananassae.fasta -l diptera_odb10 -o D.ananassae.fasta.BUSCO -m genome -c 6 -f 

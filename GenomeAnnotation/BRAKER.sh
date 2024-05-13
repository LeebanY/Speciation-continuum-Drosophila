#!/bin/bash
#SBATCH --job-name=braker    # Job name
#SBATCH --array=1-50
#SBATCH --ntasks=2                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=braker.log   # Standard output and error log
pwd; hostname; date

conda activate braker2
FILES=(/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/ASSEMBLIES/flyseq/*.fasta)
export ALIGNMENT_TOOL_PATH=/mnt/shared/scratch/lyusuf/apps/conda/envs/braker2/bin
export AUGUSTUS_BIN_PATH=/mnt/shared/scratch/lyusuf/apps/conda/envs/braker2/bin
export AUGUSTUS_SCRIPTS_PATH=/mnt/shared/scratch/lyusuf/apps/conda/envs/braker2/bin
export AUGUSTUS_CONFIG_PATH=/mnt/shared/scratch/lyusuf/apps/conda/envs/braker2/config
export GENEMARK_PATH=/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/ASSEMBLIES/BRAKER/gmes_linux_64
export CDBTOOLS_PATH=/mnt/shared/scratch/lyusuf/apps/conda/envs/braker2/bin
export DIAMOND_PATH=/mnt/shared/scratch/lyusuf/apps/conda/envs/braker2/bin
export PERL5LIB=/mnt/shared/scratch/lyusuf/apps/conda/envs/braker2/bin
export PROTHINT_PATH=/home/lyusuf/scratch/Speciation_Gimble_Drosophila_Project/ASSEMBLIES/BRAKER/ProtHint/bin
prot_evd=dmel-all-translation-r6.39.fasta

mkdir ${FILES[$SLURM_ARRAY_TASK_ID]}_braker

braker.pl --genome=${FILES[$SLURM_ARRAY_TASK_ID]} \
--prot_seq=$prot_evd \
--ALIGNMENT_TOOL_PATH=/mnt/shared/scratch/lyusuf/apps/conda/envs/braker2/bin --gff3 --verbosity=3 --cores=2 --softmasking --species=${FILES[$SLURM_ARRAY_TASK_ID]}_braker \
--workingdir=${FILES[$SLURM_ARRAY_TASK_ID]}_braker





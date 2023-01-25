#!/bin/bash

##########CONFIGURATION SLURM#####################
##nom du job##
#SBATCH --job-name=AnapSPAdes
##partition##
#SBATCH -p normal
##mail
#SBATCH --mail-user=olivier.duron@cnrs.fr
##envoi mail ##
#SBATCH --mail-type=END
##nombre de coeurs
#SBATCH -c 12
#################################################

### creation de variables
path_to_scratch='/scratch/duron_AnapSPAdes'
path_to_project='nas3:/data3/projects/GUYAVEC/'

mkdir $path_to_scratch
cd $path_to_scratch


## lancement de la boucle
reads_ech_list="BTR250 ORP110"

for reads_ech in $reads_ech_list
do

scp $path_to_project/databrut/$reads_ech-R1-trimmed.fastq.gz $path_to_scratch
scp $path_to_project/databrut/$reads_ech-R2-trimmed.fastq.gz $path_to_scratch

## assemblage SPAdes
module load bioinfo/SPAdes/3.15.3
spades.py --meta -1 $reads_ech-R1-trimmed.fastq.gz -2 $reads_ech-R2-trimmed.fastq.gz -o $reads_ech-SPAdes -t 12

rm $reads_ech-R1-trimmed.fastq.gz
rm $reads_ech-R2-trimmed.fastq.gz

scp -r $path_to_scratch/$reads_ech-SPAdes $path_to_project/Assembly_SPAdes/

done 

#!/bin/bash

##########CONFIGURATION SLURM#####################
##nom du job##
#SBATCH --job-name=AnapMEGAHIT
##partition##
#SBATCH -p normal
##mail
#SBATCH --mail-user=olivier.duron@cnrs.fr
##envoi mail ##
#SBATCH --mail-type=END
##nombre de coeurs
#SBATCH -c 8
#################################################

### creation de variables
path_to_scratch='/scratch/duron_AnapMEGAHITt'
path_to_project='nas3:/data3/projects/GUYAVEC/'

mkdir $path_to_scratch
cd $path_to_scratch

## lancement de la boucle
reads_ech_list="BTR250 OPR110"

for reads_ech in $reads_ech_list
do

scp $path_to_project/databrut/$reads_ech-R1-trimmed.fastq.gz $path_to_scratch
scp $path_to_project/databrut/$reads_ech-R2-trimmed.fastq.gz $path_to_scratch

module load bioinfo/MEGAHIT/1.2.9

megahit -1 $reads_ech-R1-trimmed.fastq.gz -2 $reads_ech-R2-trimmed.fastq.gz --k-list 21,33,41,55,77 -t 8 -o $reads_ech-metaMEGAHIT 

rm $reads_ech-R1-trimmed.fastq.gz
rm $reads_ech-R2-trimmed.fastq.gz

scp -r $path_to_scratch/$reads_ech-metaMEGAHIT  $path_to_project/Assembly_MEGAHIT/

done 


#!/bin/bash

##########CONFIGURATION SLURM#####################
##nom du job##
#SBATCH --job-name=AnaBinning
##partition##
#SBATCH -p normal
##mail
#SBATCH --mail-user=olivier.duron@ird.fr
##envoi mail ##
#SBATCH --mail-type=END
#################################################

### creation de variables
path_to_scratch='/scratch/duron_AnaBinning_BTR'
path_to_project='nas3:/data3/projects/GUYAVEC/'
path_to_shortreads='nas3:/data3/projects/GUYAVEC/databrut/'

mkdir $path_to_scratch
cd $path_to_scratch

## lancement de la boucle
reads_ech_list="BTR250"

for reads_ech in $reads_ech_list
do

scp $path_to_project/Assembly_MEGAHIT/$reads_ech-metaMEGAHIT/final.contigs.fa $path_to_scratch
scp $path_to_shortreads/$reads_ech-R1-trimmed.fastq.gz $path_to_scratch
scp $path_to_shortreads/$reads_ech-R2-trimmed.fastq.gz $path_to_scratch

## PARTIE 1
## mapper reads sur assembly.fasta
module load bioinfo/bwa/0.7.17
module load bioinfo/samtools/1.9
bwa index final.contigs.fa
bwa mem -t 4 final.contigs.fa $reads_ech-R1-trimmed.fastq.gz $reads_ech-R2-trimmed.fastq.gz | samtools sort -@ 4 -T mapped -O BAM -o $reads_ech-reads-mapped.bam
samtools index $reads_ech-reads-mapped.bam

## lancer CONCOCT
module load bioinfo/concoct/1.1.0

cut_up_fasta.py final.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

concoct_coverage_table.py contigs_10K.bed $reads_ech-reads-mapped.bam > coverage_table.tsv

concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b $reads_ech-concoctmetaMEGAHIT_output/ -t 4

merge_cutup_clustering.py $reads_ech-concoctmetaMEGAHIT_output/clustering_gt1000.csv > $reads_ech-concoctmetaMEGAHIT_output/clustering_merged.csv

mkdir $reads_ech-concoctmetaMEGAHIT_output/fasta_bins
extract_fasta_bins.py final.contigs.fa $reads_ech-concoctmetaMEGAHIT_output/clustering_merged.csv --output_path $reads_ech-concoctmetaMEGAHIT_output/fasta_bins

scp -r ./$reads_ech-concoctmetaMEGAHIT_output $path_to_project/Binning_MEGAHIT
scp $reads_ech-reads-mapped.bam $path_to_project/Binning_MEGAHIT/$reads_ech-concoctmetaMEGAHIT_output/


## PARTIE 2
#scp $path_to_project/Assembly_MEGAHIT/$reads_ech-metaMEGAHIT/final.contigs.fa $path_to_scratch
#scp $path_to_project/Binning_MEGAHIT/$reads_ech-concoctmetaMEGAHIT_output/$reads_ech-reads-mapped.bam $path_to_scratch

## faire contig.db et profile.db 
#sed '/^>/s/ .*//' final.contigs.fa > final.contigs-rename.fa
#rm final.contigs.fa  
#cp final.contigs-rename.fa final.contigs.fa

#module load system/Miniconda3/1.0
#source activate anvio-7.1

#anvi-gen-contigs-database -f final.contigs.fa -o $reads_ech-metaMEGAHIT.db --ignore-internal-stop-codons -n Binning -T 4
#anvi-run-hmms -c $reads_ech-metaMEGAHIT.db
#anvi-profile -i $reads_ech-reads-mapped.bam -c $reads_ech-metaMEGAHIT.db --min-contig-length 250 -T 4 -o $reads_ech-PROFILE --cluster-contigs

#mv $path_to_scratch/$reads_ech-PROFILE/* $path_to_scratch

## pour créer bins.txt
## The clustering_merged CSV file produced by CONCOCT needed to be exported in a tabular-delimited TEXT file. The bins' names had to be renamed to match the requirements of anvi'o :
## awk '$2="bin"$2 {print}' bins-to-format.txt > bins-renamed.txt
## The bins-renamed.txt file had to be transformed again to correspond to a tabular-delimited TEXT file, called bins.txt hereafter.

## creer collection avec les bins
#scp $path_to_project/$reads_ech-metaMEGAHIT/$reads_ech-concoctmetaMEGAHIT_output/bins.txt $path_to_scratch
#anvi-import-collection bins.txt -c $reads_ech-metaMEGAHIT.db -p PROFILE.db -C bins --contigs-mode

## attribuer taxonomy
#anvi-run-scg-taxonomy -c $reads_ech-metaMEGAHIT.db -T 2

## estimer taxonomy 
#anvi-estimate-scg-taxonomy -c $reads_ech-metaMEGAHIT.db --output-file $reads_ech-TAXONOMY.txt -p PROFILE.db -C bins --compute-scg-coverages -T 2

## summarize
#anvi-summarize -p PROFILE.db -c $reads_ech-metaMEGAHIT.db -o SUMMARY -C bins

#rm $path_to_scratch/bins.txt
#rm $reads_ech-reads-mapped.bam
#rm $reads_ech-reads-mapped.bam.bai
#rm $reads_ech-R1-trimmed.fastq.gz
#rm $reads_ech-R2-trimmed.fastq.gz
#rm final.contigs.fa 

#scp -r $path_to_scratch/* $path_to_project/Assembly_MEGAHIT/$reads_ech-metaMEGAHIT/$reads_ech-concoctmetaMEGAHIT_output/

done 

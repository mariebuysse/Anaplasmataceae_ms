# "Emerging anaplasmosis and ehrlichiosis in the Amazon biome", Buysse et al.
# Bioinformatic analyses and scripts

For this study, we aimed to retrieve the MAG (Metagenome-Assembled Genome) of three Anaplasmataceae detected in humans, sloths, and ticks. We characterized the newly obtained MAGs, then compared them to other Anaplasmataceae' genomes, especially in the search of homologs of virulence factors.
\
For each sample, a dataset of paired-end reads is available, each representing an individual metagenome. Details about the experimental and sequencing methods are available in the associated manuscript. Following analyses are based on these reads' datasets, referred as `ORP110` (reads from **X**), `BTR250` (reads from **X**), and `AcajP1` (reads from **X**).

# Step 1. Retrieving new Anaplasmataceae's MAGs
## 1.1. *De novo* assembly from the short reads' datasets: `ORP110` and `BTR250`
### 1.1.1. Trimming
```
atropos -T 4 -a file:BTR250-adaptersF -A file:BTR250-adaptersR -o BTR250-R1-trimmed.fastq.gz -p BTR250-R2-trimmed.fastq.gz -pe1 BTR250-R1.fastq.gz -pe2 BTR250-R2.fastq.gz
```
### 1.1.2. Assembly
### 1.1.2. Binning 
## 1.2. *De novo* assembly from the long reads' dataset: `AcajP1`
## 1.3. Quality check
## 1.4. Annotation 

# Step 2. MAGs' description and comparison with others genomes of Anaplasmataceae
## 2.1. Phylogenomics 
## 2.2. Gene content and nucleic sequences' similarities
### Venn diagram 
```
orthofinder -f ./FAA_Rlusi -t 5 -a 1 -S diamond ## FAA_Rlusi being a directory including .faa files of R. lusitaniae MAGs

#Ehrlichia
# On RStudio
ortho_tab <- read.table("Ehr_GeneCount.txt", sep="\t", header=TRUE, fill=TRUE) 
### list of orthologs per each organism in the dataset
Acaj <- subset(ortho_tab, ortho_tab$Acaj>0)
Acaj_list <- Acaj$Orthogroup
Ehrcani <- subset(ortho_tab, ortho_tab$Ehrcani>0)
Ehrcani_list <- Ehrcani$Orthogroup
Ehrchaff <- subset(ortho_tab, ortho_tab$Ehrchaff>0)
Ehrchaff_list <- Ehrchaff$Orthogroup
Ehrmuri <- subset(ortho_tab, ortho_tab$Ehrmuri>0)
Ehrmuri_list <- Ehrmuri$Orthogroup
Ehrrumi <- subset(ortho_tab, ortho_tab$Ehrrumi>0)
Ehrrumi_list <- Ehrrumi$Orthogroup
Orthologs <- list(Acaj = Acaj_list, Ehrcani = Ehrcani_list, Ehrchaff = Ehrchaff_list, Ehrmuri = Ehrmuri_list, Ehrrumi = Ehrrumi_list) ## create a list of lists

library(VennDiagram)
set.seed(1)
venn.diagram(Orthologs, filename="Ehr_comparison.png", imagetype = "png", height=2000, width=2000, cex=0.8, cat.cex=0.8, fill=c("#22780F", "#CECECE", "#798081","#C1BFB1", "#BBACAC"), col="black", lwd=1, cat.dist=0.25) # change height, width, color 

#Anaplasma
# On RStudio
ortho_tab <- read.table("Anap_GeneCount.txt", sep="\t", header=TRUE, fill=TRUE) 
### list of orthologs per each organism in the dataset
Amargi <- subset(ortho_tab, ortho_tab$Amargi>0)
Amargi_list <- Amargi$Orthogroup
Aphago <- subset(ortho_tab, ortho_tab$Aphago>0)
Aphago_list <- Aphago$Orthogroup
Aplaty <- subset(ortho_tab, ortho_tab$Aplaty>0)
Aplaty_list <- Aplaty$Orthogroup
BTR250 <- subset(ortho_tab, ortho_tab$BTR250>0)
BTR250_list <- BTR250$Orthogroup
ORP110 <- subset(ortho_tab, ortho_tab$ORP110>0)
ORP110_list <- ORP110$Orthogroup
Orthologs <- list(ORP110 = ORP110_list, Amargi = Amargi_list, Aphago = Aphago_list, Aplaty = Aplaty_list, BTR250 = BTR250_list) ## create a list of lists

library(VennDiagram)
set.seed(1)
venn.diagram(Orthologs, filename="Anap_comparison.png", imagetype = "png", height=2000, width=2000, cex=0.8, cat.cex=0.8, fill=c("#CF0A1D", "#CECECE", "#798081","#C1BFB1","#318CE7"), col="black", lwd=1, cat.dist=0.25) # change height, width, color 
```

### ANI
With pyani's script `average_nucleotide_identity.py` (https://github.com/widdowquinn/pyani/blob/master/README_v_0_2_x.md)
```
average_nucleotide_identity.py -i Genomes_fastANI/ -o pyani_results -g
```

## 2.3. Detection of virulence factors 

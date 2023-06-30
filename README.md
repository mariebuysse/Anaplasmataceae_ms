# Anaplasma_analyses

#### Raw data
```
nas3:/data3/projects/GUYAVEC/databrut
```
Echantillons : BTR250 & ORP110
\
Run : flow cell SP en paired-end 150 nt

# Quality check
With FASTQC, with control of adapters:
```
fastqc BTR250-R1.fastq.gz
```
If trimming needed: 
```
atropos -T 4 -a file:BTR250-adaptersF -A file:BTR250-adaptersR -o BTR250-R1-trimmed.fastq.gz -p BTR250-R2-trimmed.fastq.gz -pe1 BTR250-R1.fastq.gz -pe2 BTR250-R2.fastq.gz
```
If need to decompress:
```
gzip -d fichier.gz
tar -xvf fichier.tar
unzip fichier.zip
tar xjf fichier.tar.bz2
tar xvfz fichier.tar.gz
tar xf fichier.tar.xz
```
Check with FASTQC again:
```
fastqc BTR250-R1-trimmed.fastq.gz
```
Check the optimal kmer size for assembly:
```
kmergenie -t 6 BTR250-R1-trimmed.fastq.gz
```

# Assembly
See details in scripts for SPAdes and MEGAHIT

# Check quality assembly
With QUAST:
```
quast.py *.fa -o ../2.Quality_assembly/QUAST_RESULTS_metaMEGAHIT
```
```
Assembly                    BRT250-final.contigs  ORP110-final.contigs
# contigs (>= 0 bp)         1636409               1168235             
# contigs (>= 1000 bp)      386703                488189              
# contigs (>= 5000 bp)      150203                165797              
# contigs (>= 10000 bp)     64913                 61294               
# contigs (>= 25000 bp)     7886                  4858                
# contigs (>= 50000 bp)     448                   121                 
Total length (>= 0 bp)      2720728327            2752617663          
Total length (>= 1000 bp)   2290459576            2494217607          
Total length (>= 5000 bp)   1714251326            1702229590          
Total length (>= 10000 bp)  1108439149            968353092           
Total length (>= 25000 bp)  263096035             153286388           
Total length (>= 50000 bp)  26922741              6975055             
# contigs                   561790                631882              
Largest contig              128915                102827              
Total length                2411301713            2596898192          
GC (%)                      39.91                 40.76               
N50                         9078                  7456                
N75                         4320                  3694                
L50                         75123                 99649               
L75                         170458                222791              
# N's per 100 kbp           0.00                  0.00
```
Results: let's try without reads from Eukaryota (remove with Kraken) because of the low size of largest contig

With Bandage:
```
module load bioinfo/MEGAHIT/1.2.9
megahit_toolkit contig2fastg 21 k21.contigs.fa > BTR250_k21.fastg
megahit_toolkit contig2fastg 77 k77.contigs.fa > BTR250_k77.fastg
```

# Binning
Test binning on metaMEGAHIT on BTR and ORP final.contigs.fa with CONCOCT coupled with anvi'o taxonomy tools - in interactive mode (see scripts)

The `clustering_merged` CSV file produced by CONCOCT needed to be exported in a tabular-delimited TEXT file. The bins' names had to be renamed to match the requirements of anvi'o :
```
awk '$2="bin"$2 {print}' bins-to-format.txt > bins-renamed.txt
```
The `bins-renamed.txt` file had to be transformed again to correspond to a tabular-delimited TEXT file, called `bins.txt` hereafter.

Results of comparison meta and filtered: meta then binning CONCOCT-anvio

```
Assembly                    BTR-bin23ANA-contigs  ORP-bin116ANA-contigs
# contigs (>= 0 bp)         81                    117                  
# contigs (>= 1000 bp)      81                    117                  
# contigs (>= 5000 bp)      47                    62                   
# contigs (>= 10000 bp)     32                    42                   
# contigs (>= 25000 bp)     16                    12                   
# contigs (>= 50000 bp)     5                     0                    
Total length (>= 0 bp)      1176232               1187243              
Total length (>= 1000 bp)   1176232               1187243              
Total length (>= 5000 bp)   1122570               1068187              
Total length (>= 10000 bp)  1012570               917159               
Total length (>= 25000 bp)  730933                441829               
Total length (>= 50000 bp)  364916                0                    
# contigs                   81                    117                  
Largest contig              112391                49298                
Total length                1176232               1187243              
GC (%)                      50.62                 49.25                
N50                         30968                 20180                
N75                         18795                 11467                
L50                         11                    19                   
L75                         24                    40                   
# N's per 100 kbp           0.00                  0.00
```

# Venn diagram 
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

# ANI
With pyani's script `average_nucleotide_identity.py` (https://github.com/widdowquinn/pyani/blob/master/README_v_0_2_x.md)
```
average_nucleotide_identity.py -i Genomes_fastANI/ -o pyani_results -g
```

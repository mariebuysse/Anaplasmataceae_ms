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

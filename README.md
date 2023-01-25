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

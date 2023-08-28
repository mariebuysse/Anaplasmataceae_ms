# "Emerging anaplasmosis and ehrlichiosis in the Amazon biome", Buysse et al.
# Bioinformatic analyses and scripts

For this study, we aimed to retrieve the MAG (Metagenome-Assembled Genome) of three Anaplasmataceae detected in humans, sloths, and ticks in French Guiana. We characterized the newly obtained MAGs, then compared them to other representative *Anaplasma* and *Ehrlichia* genomes, especially in the search of homologs of virulence effectors.
\
In details, we sequenced metagenomes from three infected samples:
1. A human blood sample positive for *Candidatus* Anaplasma sparouinense, hereafter referred as `Sparouine`
2. A blood sample of a brown-throated three-toed sloth *Bradypus variegatus* positive for *Ca.* Anaplasma amazonensis, hereafter referred as `PetitSaut`
3. An *Amblyomma cajennense* tick sample positive for an Ehrlichia sp., hereafter designated as the putative new species *Ca.* Ehrlichia cajennense and referred as `Matoury`.
\
For `Aspa-Sparouine` and `Aama-PetitSaut`, we obtained paired-end short reads' datasets and for `Ecaj-Matoury` we obtained a long reads' dataset (see Material and method section within the manuscript for details). 
\
Hereafter,
```
ech="Sparouine PetitSaut Matoury"
```

# Step 1. Retrieving new Anaplasmataceae's MAGs
## 1.1. *De novo* assembly from the short reads' datasets: `Sparouine` and `PetitSaut`
### 1.1.1. Trimming
```
atropos -T 4 -a file:$ech-adaptersF -A file:$ech-adaptersR -o $ech-R1-trimmed.fastq.gz -p $ech-R2-trimmed.fastq.gz -pe1 $ech-R1.fastq.gz -pe2 $ech-R2.fastq.gz
```

### 1.1.2. Assembly
Reads were assembled using MEGAHIT (https://github.com/voutcn/megahit, Li D., Liu C-M., Luo R., Sadakane K., and Lam T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics. doi: 10.1093/bioinformatics/btv033):
```
megahit -1 $ech-R1-trimmed.fastq.gz -2 $ech-R2-trimmed.fastq.gz --k-list 21,59,77,99 -t 6 -o $ech-metaMEGAHIT
```

### 1.1.3. Binning 
*Anaplasma* MAGs were retrieved from assemblies using `CONCOCT` (https://github.com/BinPro/CONCOCT, Alneberg J., Smári Bjarnason B., de Bruijn I., Schirmer M., Quick J., Ijaz U.Z., Lahti L., Loman N.J., Andersson A.F., & Quince C., (2014) Binning metagenomic contigs by coverage and composition, Nature Methods. doi: 10.1038/nmeth.3103) and the `anvi'o` pipeline (https://anvio.org/, Eren A.M., Kiefl E., Shaiber A. et al., (2021) Community-led, integrated, reproducible multi-omics with anvi’o. Nature Microbiology. doi: 10.1038/s41564-020-00834-3). 
First, the contigs were renamed to match the requirements of `anvi'o`. 
```
sed '/^>/s/ .*//' $ech-final.contigs.fa > $ech-final.contigs-rename.fa
rm $ech-final.contigs.fa
mv $ech-final.contigs-rename.fa $ech-final.contigs.fa
```
The contigs were binned using `CONCOCT`: 
```
bwa index $ech-final.contigs.fa
bwa mem -t 4 $ech-final.contigs.fa $ech-R1.fastq.gz $ech-R2.fastq.gz | samtools sort -@ 4 -T mapped -O BAM -o $ech-reads-mapped.bam
samtools index $ech-reads-mapped.bam

cut_up_fasta.py $ech-final.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed $ech-reads-mapped.bam > coverage_table.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b $ech-concoctMEGAHIT_output/ -t 4
merge_cutup_clustering.py $ech-concoctMEGAHIT_output/clustering_gt1000.csv > $ech-concoctMEGAHIT_output/clustering_merged.csv
mkdir $ech-concoctMEGAHIT_output/fasta_bins
extract_fasta_bins.py $ech-final.contigs.fa $ech-concoctMEGAHIT_output/clustering_merged.csv --output_path $ech-concoctMEGAHIT_output/fasta_bins
```
The `clustering_merged` CSV file produced by CONCOCT needed to be exported in a tabular-delimited TEXT file. 
The bins' names had to be renamed to match the requirements of `anvi'o` : 
```
awk '$2="bin"$2 {print}' bins-to-format.txt > bins-renamed.txt
```
The `bins-renamed.txt` file had to be transformed again to correspond to a tabular-delimited TEXT file, called `bins.txt` hereafter.

The bins were taxonomically assigned using the `anvi'o` pipeline below: 
```
## To create the contigs database
anvi-gen-contigs-database -f $ech-final.contigs.fa -o $ech-metaMEGAHIT.db --ignore-internal-stop-codons -n Binning -T 4
anvi-run-hmms -c $ech-metaMEGAHIT.db

## To create the profile database
anvi-profile -i $ech-reads-mapped.bam -c $ech-metaMEGAHIT.db --min-contig-length 250 -T 4 -o $ech-PROFILE --cluster-contigs

## To import the bins as a collection
anvi-import-collection bins.txt -c $ech-metaMEGAHIT.db -p $ech-PROFILE/PROFILE.db -C bins --contigs-mode

## To assign the bins and visualize the results
anvi-run-scg-taxonomy -c $ech-metaMEGAHIT.db -T 2
anvi-estimate-scg-taxonomy -c $ech-metaMEGAHIT.db --output-file $ech-TAXONOMY.txt -p $ech-PROFILE/PROFILE.db -C bins --compute-scg-coverages -T 2
anvi-summarize -p $ech-PROFILE/PROFILE.db -c $ech-metaMEGAHIT.db -C bins -o $ech-SUMMARY 
```
`$ech-SUMMARY` file (html format) was manually checked to identify the bin number assigned as *Anaplasma spp.*. No redundacy level required to refine the raw bins. 
**After this step, the bins were considered as *Ca.* Anaplasma MAGs, designated as `Aspa-Sparouine` and `Aama-PetiSaut`.**

## 1.2. *De novo* assembly from the long reads' dataset: `Matoury`
*De novo* assembly of `Matoury` dataset was performed from long-reads (MinION, Oxford Nanopore) using `Flye` (https://github.com/fenderglass/Flye, M.K., D.M.B., B.B., A.G., M.R., S.B.S., K.K., J.Y., E.P., T.P.L.S. and P.A.P., (2020) metaFlye: scalable long-read metagenome assembly using repeat graphs, Nature Methods. doi:s41592-020-00971-x), as follows:
```
module load bioinfo/Flye/2.4.1
flye --nano-raw $ech_porechopped_all.fq.gz --out-dir ./$ech-Flye --threads 12 --iterations 5 --meta --min-overlap 8000 --debug --genome-size 200000
```

The resulting file `assembly.fasta` is then corrected using only the long-reads' dataset using `Medaka` Oxford Nanopore tool (https://github.com/nanoporetech/medaka):
```
mkdir ./Flye-Medaka
samtools faidx assembly.fasta
module load bioinfo/medaka/1.5
medaka_consensus -i $ech_porechopped_all.fq.gz -d assembly.fasta -m r941_min_fast_g303 -o Flye-Medaka -t 4
```

The identification of the *Ca.* Ehrlichia cajennense MAG is based on the taxonomic assignation of the 16S rDNA sequence of a circular contig visualized with `Bandage` (https://rrwick.github.io/Bandage/, Wick R.R, Schultz M.B., Zobel J., Holt K.E. (2015) Bandage: Interactive visualization of de novo genome assemblies. Bioinformatics. doi: 10.1093/bioinformatics/btv383) (assignation based on the online NCBI BLAST tool). 

**After this step, the contig is considered as *Ca.* Ehrlichia cajennense MAG, designated as `Ecaj-Matoury`.**


## 1.3. Quality check
Quality and multiple statistics were accessed using `miComplete` (https://pypi.org/project/micomplete/, Hugoson E., Lam W.T., Guy L. (2020) miComplete: weighted quality evaluation of assembled microbial genomes. Bioinformatics. doi: 10.1093/bioinformatics/btz664) and `Quast` (https://github.com/ablab/quast, Gurevich A., Saveliev V., Vyahhi N., Tesler G. (2013) QUAST: quality assessment tool for genome assemblies. Bioinformatics. doi: 10.1093/bioinformatics/btt086).
```
## With Quast
quast.py ./MAGs/* -o QUAST #a directory containing all MAGs 
## With miComplete
miComplete set.tab --hmms Bact105 #set.tab a tabular separated file containing per line both a path to each MAG and the type (here fna)
```
```
Results are:

```

## 1.4. Annotation 
The MAGs were annotated using `Prokka` (https://github.com/tseemann/prokka, Seemann T. (2014) Prokka: rapid prokaryotic genome annotation. Bioinformatics. doi: 10.1093/bioinformatics/btu153).
```
prokka $ech-MAG.fasta --locustag $ech --prefix $ech --outdir PROKKA-$ech --rfam --compliant
````

## 1.5. Genome visualization
MAG representation was performed using `CGview` (https://github.com/paulstothard/cgview, Stothard P., Wishart D.S. (2005) Circular genome visualization and exploration using CGView. Bioinformatics. doi: 10.1093/bioinformatics/bti054):
```
perl ~/Tools_starters/cgview_xml_builder.pl -sequence ./$ech.gbf -output $ech.xml -gc_skew F -gc_content F -size large-v2 -gene_decoration arc
java -jar ~/Tools_starters/cgview.jar -i $ech.xml -o map_$ech.png -f png
```


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

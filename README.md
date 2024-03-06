# Novel Anaplasma and Ehrlichia bacteria in humans, wildlife, and ticks in the Amazon rainforest, Buysse et al.
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
## With Quast
Assembly                    Ecaj-Matoury  Aama-PetiSaut   Aspa-Sparouine 
# contigs (>= 0 bp)         1        81       117    
# contigs (>= 1000 bp)      1        81       117    
# contigs (>= 5000 bp)      1        47       62     
# contigs (>= 10000 bp)     1        32       42     
# contigs (>= 25000 bp)     1        16       12     
# contigs (>= 50000 bp)     1        5        0      
Total length (>= 0 bp)      1177323  1176232  1187243
Total length (>= 1000 bp)   1177323  1176232  1187243
Total length (>= 5000 bp)   1177323  1122570  1068187
Total length (>= 10000 bp)  1177323  1012570  917159 
Total length (>= 25000 bp)  1177323  730933   441829 
Total length (>= 50000 bp)  1177323  364916   0      
# contigs                   1        81       117    
Largest contig              1177323  112391   49298  
Total length                1177323  1176232  1187243
GC (%)                      32.08    50.62    49.25  
N50                         1177323  30968    20180  
N75                         1177323  18795    11467  
L50                         1        11       19     
L75                         1        24       40     
# N's per 100 kbp           0.00     0.00     0.00   

## With miComplete
Name	Length	GC-content	Present Markers	Completeness	Redundancy	Contigs	N50	L50	N90	L90
AcajEhr	1177323	32.08	95	0.9048	1.0842	1	1177323	1	117732311410	
BTR250	1176232	50.61	100	0.9524	1.0000	81	30968	11	7853	38
ORP110	1187243	49.23	101	0.9619	1.0099	117	20180	19	4879	63
```

## 1.4. Annotation and pseudogene identification
The MAGs were annotated using `Prokka` (https://github.com/tseemann/prokka, Seemann T. (2014) Prokka: rapid prokaryotic genome annotation. Bioinformatics. doi: 10.1093/bioinformatics/btu153), then outpus are analyzed through `pseudofinder` (https://github.com/filip-husnik/pseudofinder, Syberg-Olsen M.J., Graber A.I., Keeling P.J., McCutcheon J.P., Husnik F. (2022) Pseudofinder: detection of pseudogenes in prokaryotic genomes. Molecular Biology and Evolution. doi: 10.1093/molbev/msac153.). Previously, database based on `$ech.fna` of each MAG were created using `diamond` (https://github.com/bbuchfink/diamond, Buchfink B., Reuter K., Drost H.G. (2021) Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods. doi:10.1038/s41592-021-01101-x). 
```
prokka $ech-MAG.fasta --locustag $ech --prefix $ech --outdir PROKKA-$ech --rfam --compliant
diamond makedb --in $ech.fna -d $ech_db
python3 /pseudofinder-1.0/pseudofinder.py annotate --genome $ech.gbk --outprefix $ech --diamond -db $ech_db
```

## 1.5. Genome visualization
MAG representation was performed using `CGview` (https://github.com/paulstothard/cgview, Stothard P., Wishart D.S. (2005) Circular genome visualization and exploration using CGView. Bioinformatics. doi: 10.1093/bioinformatics/bti054):
```
perl ~/Tools_starters/cgview_xml_builder.pl -sequence ./$ech.gbf -output $ech.xml -gc_skew F -gc_content F -size large-v2 -gene_decoration arc
java -jar ~/Tools_starters/cgview.jar -i $ech.xml -o map_$ech.png -f png
```


# Step 2. MAGs' description and comparison with others genomes of Anaplasmataceae
## 2.1. Phylogenomics 
First, single-copy orthologs (SCO) were identified using OrthoFinder (https://github.com/davidemms/OrthoFinder, Emms D.M. and Kelly S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biology. doi: 10.1186/s13059-019-1832-y) from a set of specimens' genomes chosen to study the phylogenetic relationships between the obtained Spiroplasma MAGs and other Spiroplasma representatives:
```
orthofinder -f ./OrthoFinder_genomes/ -t 4 -S blast ## OrthoFinder_genomes being a directory including all .faa files of specimens of interest
```
For each SCO, sequences were individually aligned using mafft (https://github.com/GSLBiotech/mafft, Katoh K. and Standley D.M. (213) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution. doi: 10.1093/molbev/mst010):
```
for file in /Single_Copy_Orthologue_Sequences/*
do mafft "$file" > "$file"
done
```
For each SCO, ambigious hypervariable regions were removed using trimAl (https://github.com/inab/trimal, Capella-Gutiérrez S., Silla-Martínez J.M., and Gabaldón T. (2009) trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics. doi: 10.1093/bioinformatics/btp348):
```
cp ./Single_Copy_Orthologue_Sequences/*_align.fasta ./Single_Copy_Orthologue_Sequences_trimal/
for file in /Single_Copy_Orthologue_Sequences_trimal/*
do trimal -in "$file" -out "$file" -fasta -gt 1 -cons 50
done
```
Then, all SCO sequences were concatenated using Amas (https://github.com/marekborowiec/AMAS, Borowiec M.L. (2016) AMAS: a fast tool for alignment manipulation and computing of summary statistics. PeerJ. doi: 10.7717/peerj.1660) in a single file:
```
for file in /Single_Copy_Orthologue_Sequences_trimal/*
do awk '/^>/{print ">organism" ++i; next}{print}' < "$file" > "${file%_align.fasta}_rename.fasta"
done
cp ./Single_Copy_Orthologue_Sequences_trimal/*_rename.fasta ./Single_Copy_Orthologue_Sequences_AMAS/
AMAS.py concat -f fasta -d aa --in-files ./Single_Copy_Orthologue_Sequences_AMAS/*.fasta
```
Substitution models were evaluated using modeltest-ng (https://github.com/ddarriba/modeltest, Darriba D., Posada D., Kozlov A.M., Stamatakis A., Morel B., and Flouri T. (2020) ModelTest-NG: A new and scalable tool for the selection of DNA and protein evolutionary models. Molecular Biology and Evolution. doi: 10.1093/molbev/msz189) in order to determinate the appropriate substitution model (according to AICc criterion) to use for the phylogenetic tree construction with RAxML (https://github.com/stamatak/standard-RAxML, Stamatakis A. (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics. doi: 10.1093/bioinformatics/btu033):
```
modeltest-ng -i SCO_concatenated.faa -p 12 -T raxml -d aa
raxmlHPC-PTHREADS -T 8 -f a -s SCO_concatenated.faa -n phylo -m PROTGAMMAIJTT -x 1234 -# 1000 -p 1234
```
The phylogenetic tree was visualized and modified using figtree (<https://github.com/rambaut/figtree/>).

## 2.2. Gene content and nucleic sequences' similarities
### Venn diagram 
Genome content (orthologs) between MAGs and representatives' genomes of each genus were compared using `OrthoFinder` and the following script to produce a Venn Diagramm on R (v3.6.3) (https://www.R-project.org/, R Core Team. (2020) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria) and RStudio (http://www.rstudio.com/, RStudio Team. (2020) RStudio: Integrated Development for R. RStudio, PBC, Boston, MA):
```
orthofinder -f ./FAA_files -t 5 -a 1 -S diamond ## FAA_files being a directory including .faa files of all MAg and genome to compare

#Ehrlichia
# On RStudio
ortho_tab <- read.table("Ehr_GeneCount.txt", sep="\t", header=TRUE, fill=TRUE) 
### list of orthologs per each organism in the dataset
Ecaj-Matoury <- subset(ortho_tab, ortho_tab$Acaj>0)
Ecaj-Matoury_list <- Ecaj-Matoury$Orthogroup
Ehrcani <- subset(ortho_tab, ortho_tab$Ehrcani>0)
Ehrcani_list <- Ehrcani$Orthogroup
Ehrchaff <- subset(ortho_tab, ortho_tab$Ehrchaff>0)
Ehrchaff_list <- Ehrchaff$Orthogroup
Ehrmuri <- subset(ortho_tab, ortho_tab$Ehrmuri>0)
Ehrmuri_list <- Ehrmuri$Orthogroup
Ehrrumi <- subset(ortho_tab, ortho_tab$Ehrrumi>0)
Ehrrumi_list <- Ehrrumi$Orthogroup
Orthologs <- list(Ecaj-Matoury = Ecaj-Matoury_list, Ehrcani = Ehrcani_list, Ehrchaff = Ehrchaff_list, Ehrmuri = Ehrmuri_list, Ehrrumi = Ehrrumi_list) ## create a list of lists
library(VennDiagram)
set.seed(1)
venn.diagram(Orthologs, filename="Ehr_comparison.png", imagetype = "png", height=2000, width=2000, cex=0.8, cat.cex=0.8, col="black", lwd=1, cat.dist=0.25)

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
Aama-PetitSaut <- subset(ortho_tab, ortho_tab$Aama-PetitSaut>0)
Aama-PetitSaut_list <- Aama-PetitSaut$Orthogroup
Aspa-Sparouine <- subset(ortho_tab, ortho_tab$Aspa-Sparouine>0)
Aspa-Sparouine_list <- Aspa-Sparouine$Orthogroup
Orthologs <- list(Aspa-Sparouine = Aspa-Sparouine_list, Amargi = Amargi_list, Aphago = Aphago_list, Aplaty = Aplaty_list, Aama-PetitSaut = Aama-PetitSaut_list) ## create a list of lists
library(VennDiagram)
set.seed(1)
venn.diagram(Orthologs, filename="Anap_comparison.png", imagetype = "png", height=2000, width=2000, cex=0.8, cat.cex=0.8, col="black", lwd=1, cat.dist=0.25)
```

### ANI
Nucleic sequences' similarity matrix was calculated between MAGs and representatives' genomes of each genus using `pyani`'s script `average_nucleotide_identity.py` (https://github.com/widdowquinn/pyani/blob/master/README_v_0_2_x.md) as follows:
```
average_nucleotide_identity.py -i Genomes_fastANI/ -o pyani_results -g ## Genomes_fastANI a directory containing MAGs and genomes of Anaplasmataceae to compare (.fna)
```

## 2.3. Detection of virulence factors 
The presence and completeness of genes involved in pathways of interest were investigated by two complementary tools to increase detection of both (potential) functionnal and pseudogenized genes:
1. Using OrthoFinder
2. Using BLASTn, BLASTp, and tBLASTn.
```
## method 1
orthofinder -f ./FAA_files/ -t 4 -S blast ## FAA_files being a directory with .faa files from all species of interest, including either the files "Multiquery_Anaplasma-virulence_prot.fasta" or "Multiquery_Ehrlichia-virulence_prot.fasta"

## method 2
### BLASTn
makeblastdb -in $ech-genome.fasta -dbtype nucl -out $ech_db
blastn -query Multiquery_Anaplasma-virulence_gene.fasta -outfmt "6 qseqid sseqid sseq qlen pident nident mismatch evalue sstart send gapopen" -perc_identity 30 -out queries-gene-blastn_vs_$ech.out -db $ech_db -num_threads 6 ## same with Multiquery_Ehrlichia-virulence_gene.fasta

### BLASTp
makeblastdb -in $ech-genome.faa -dbtype prot -out $ech_db
blastp -query Multiquery_Anaplasma-virulence_prot.fasta -outfmt "6 qseqid sseqid sseq qlen pident nident mismatch evalue sstart send gapopen" -evalue 1e-10 -out queries-protein-blastp_vs_$ech.out -db $ech_db -num_threads 6 ## same with Multiquery_Ehrlichia-virulence_prot.fasta

### tBLASTn
makeblastdb -in $ech-genome.fasta -dbtype nucl -out $ech_db
tblastn -query Multiquery_Anaplasma-virulence_prot.fasta -outfmt "6 qseqid sseqid sseq qlen pident nident mismatch evalue sstart send gapopen" -evalue 1e-10 -out queries-prot-tblastn_vs_$ech.out -db $ech_db -num_threads 6 ## same with Multiquery_Ehrlichia-virulence_prot.fasta
```

The visualization of the structure of Anaplasma msp2 operon and Ehrlichia omp-1/p28 proteins were produced using genoplotR.

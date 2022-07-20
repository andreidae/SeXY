# SeXY

SeXY is a simple sex-identification method for taxa lacking assembled conspecific sex chromosomes, applicaple to shotgun data from mammal species with a heterogametic and a homogametic sex and potentially to other species in which the target and reference species share the same sex determination system.

It uses a synteny-based approach to identify putative X-linked scaffolds in a given reference assembly, and take advantage of the expectation that males have half the amount of X-chromosome genetic material compared to females.


## The seXY method requires: 
1. Raw sequencing reads of the target individual, 
2. An assembled genome of the target species/or closely related species or `reference genome assembly`. 
3. Assembled X and Y chromosomes from a related chromosome-level assembly or `reference sex-chromosome assembly`.


## In general, the steps of SeXY are:

A. Use the reference sex-chromosome assemblies (X and Y) to identify putative sex-linked scaffolds in the reference genome assembly via synteny.

B. Extract scaffolds aligning to the X reference sex chromosome assembly and remove pseudoautosomal regions, i.e. short homology regions between the X and Y chromosome. 

C. Extract autosomal scaffolds, i.e. scaffolds not aligning to either the X or Y reference sex-chromosome assembly. 

D. Map raw shotgun reads to the reference genome assembly.

E. Calculate the mean coverage of ten million sites randomly sampled across the extracted X-linked scaffolds. 

F. Calculate the mean coverage of ten million sites randomly sampled across the extracted autosomal scaffold.

G. Calculate the ratio of X-linked scaffold mean coverage to Autosomal scaffold mean coverage (X:A ratio).

H. Repeat coverage calculations (E-G) ten times. 

I. Calculate the mean and standard deviations of the X:A ratio.


## Pipeline steps

### Required software
Satsuma synteny (http://satsuma.sourceforge.net/)

BEDtools (https://bedtools.readthedocs.io/en/latest/)

SAMtools (http://samtools.sourceforge.net/)


### Additional software
BBmap (https://jgi.doe.gov/data-and-tools/bbtools/)


### Download the reference genome assembly
See frequently asked questions


### Download the reference sex-chromosome assembly
See frequently asked questions


### A. Run satsuma synteny to find which scaffolds align with sex chromosomes.
Before running synteny, concatenate the RefX and RefY assemblies in one file using the command:

`cat RefX.fasta RefY.fasta > RefX_RefY.fasta`

The RefX and RefY can also be kept separate, but you then need to run two synteny analyses.


Use this command to run satsuma synteny:  

`sh SatsumaSynteny.sh $1 $2 $3 $4 $5`

```
Replace each number with the following variables
$1 - Threads
$2 - Query sequence (reference genome assembly, RefGen)
$3 - Target sequence (reference sex-chromosome assembly, RefX and RefY)
$4 - Output_directory
$5 - Satsuma directory

```
 - To reduce memory and time requirements you can remove all scaffolds <10kb from the reference genomes prior to alignment. This can be done easily using reformat.sh from the BBmap toolsuite (Bushnell, 2014), for example:

`reformat.sh in=RefGen.fasta out=RefGen_10kb.fasta minlength=10000`

where RefGen.fasta is your original RefGen and RefGen_10kb.fasta is your RefGen after removing scaffolds <10 kb

- To run reformat.sh you can copy the file to your working directory from BBmap toolsuit GitHub
https://github.com/BioInfoTools/BBMap/blob/master/sh/reformat.sh



### B-C Create bed files of X-scaffolds and autosomal scaffolds
- Extract regions mapping to X and Y from the satsuma output independently

`grep chromosome_Y satsuma_summary.chained.out | awk '{print $4"\t"$5"\t"$6}' > Y.bed`

`grep chromosome_X satsuma_summary.chained.out | awk '{print $4"\t"$5"\t"$6}' > X.bed`

Note: the terms chromosome_Y and chromosome_X may vary and will depend on the fasta file header of your reference sex-chromosome assembly.

- Remove overlapping regions from the X bed file (putatively pseudoautosomal regions)

`grep -v -f <(bedtools intersect -a X.bed -b Y.bed) X.bed > X_trim.bed`

- Extract scaffolds not mapping to X and Y from the satsuma output independently
 
`samtools faidx RefGEN_10kb.fasta`

`grep -v -f <(cat Y.bed X.bed | cut -f 1 | sort | uniq) RefGEN_10kb.fai | awk '{print $1"\t1\t"$2}'  > Autosomes.bed`


### D. Map raw shotgun reads to the reference genome assembly
Map your raw shotgun reads to RefGen. The RefGen includes both autosome- and sex- chromosome scaffolds. If you remove the short scaffolds <10 kb, use RefGen_10kb.fasta. The output from the mapping is a bam file. See frequently asked questions.


### E.I. Calculate depths of bam files

-  Randomly sample 1M sites 10x and calculate depth from said sites

`sh Coverage_calculation.sh Autosomes.bed X_trim.bed Output_directory Bamfile Output_prefix`

This will result in a txt file called `Output_directory/Outputprefix_ratios.txt` where the 10 X:A ratios can be found


##  Optional step
### Downsampling the bam files
#### (Step between D. “Map raw reads to the RefGen” and E. “Calculate depths of bam files”)
This is an optional step. It was used for testing how the number of reads in a bam file can influence sex identification. 

It requires BBmap toolsuite (Bushnell, 2014) and samtools.
We show an example for one of the species and RefGen analysed in our study

`cat commands_bbmap_dog.txt` 

``` 
in commands_bbmap_dog.txt
./bbMAP_PB_DOG.sh sampleID_1
./bbMAP_PB_DOG.sh sampleID_2
./bbMAP_PB_DOG.sh sampleID_3
./bbMAP_PB_DOG.sh sampleID_4
./bbMAP_PB_DOG.sh sampleID_5
./bbMAP_PB_DOG.sh sampleID_6
./bbMAP_PB_DOG.sh sampleID_7
./bbMAP_PB_DOG.sh sampleID_8
./bbMAP_PB_DOG.sh sampleID_9
./bbMAP_PB_DOG.sh sampleID_10

```

```
in bbMAP_PB_DOG.sh
src=$WDIR/ref_DOG_10kb
ref=DOG

for i in 100000 50000 10000 5000 2500 1000; do for j in 1 2 3 4 5;

do

reformat.sh \
	in=${src}/${1}.${ref}_10kb.bam \
	out=${1}.${ref}_10kb_${i}_${j}.bam \
	samplereadstarget=${i} \
	sampleseed=1${j}0

  done; done

```

The number of reads in the bam file can be checked with:

` samtools view -c .bam`


## Frequently asked questions

**1. How to download the reference genome assembly?**


You can download a genome from a genome assembly depository such as NCBI or DNAzoo. The steps for NCBI are:
- Go to NCBI: https://www.ncbi.nlm.nih.gov/
- In the search option select: “Assembly” and write the species name.
- Select the assembly species and version of interest.
- On the right panel, “Access the data”, click on “FTP directory for RefSeq assembly”.
- Copy the link address of the genome in .fna.gz format.
- You can use this address to download the genome using the command line `wget [fna.gz genome link]`. 

E.g. Beluga.v3 assembly reference genome: 

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/288/925/GCF_002288925.2_ASM228892v3/GCF_002288925.2_ASM228892v3_genomic.fna.gz`


**2. How to obtain a reference sex-chromosome assembly?**

This can either be done by downloading the sex chromosomes directly or by downloading a chromosome level assembly and extracting the sex chromosomes.
To download the sex chromosomes directly, follow the same steps as to with the reference genome. In “FTP directory for RefSeq assembly” there is usually a folder called assembly_structure in which the individual chromosomes can be downloaded in fasta format `wget [fna.gz sex-chromosome link]`. 

E.g. Xchr (RefX) from the cow:

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chrX.fna.gz`

To extract the sex chromosomes from RefGen one first needs to find the fasta header for the chromosomes, once that is available samtools faidx can be used. E.g. 

`samtools faidx reference.fasta Xheader > Xchromosome.fasta`


**3. How to map raw reads to the reference genome assembly?**

You can use a developed pipeline such as PALEOMIX (https://paleomix.readthedocs.io/en/stable/). This pipeline is designed to process High Throughput Sequencing data. Was originally designed with ancient DNA data but can also be used for processing modern samples. 

**4. How do I know that the mitochondrial genome is not included in the reference genome assembly?**

Most reference genome assemblies exclude the mitochondrial genome. You should be able to see that information in the description of the assembly in e.g. NCBI. 
If it is uncertain whether the mitochondrial is included or not, you can map your reference genome to the mitochondrial genome and remove those reads, or you can exclude all scaffolds <18 kb. E.g.

`reformat.sh in=RefGen.fasta out=RefGen_18kb.fasta minlength=18000`

**5. Does SeXY work on other species besides mammals?**

Although we assessed the method based on XY sex chromosomes (as in mammals), the method can in theory be applied to any species with a heterogametic and a homogametic sex (e.g, birds, and some reptiles, fish, and insects), and in which the target and reference species share the same sex determination system. However, it will requiere additional evaluations.
Here we show the results on an avian species (snow goose) that we evalueted.

Table 1: Results of the seXY pipeline in determining the sex of an avian species

|   |  |  |
| ------------- | ------------- | ------------- |
| **Species** | *Anser caerulescens* | *Anser caerulescens* | 
| **Common name** | Snow goose | Snow goose |
| **Expected sex** | M | F |
| **RefZW** | Chicken | Chicken |
| **RefGen** | *Anser brachyrhynchus* | *Anser brachyrhynchus* |
| **RefGen N50 (bp)** | 4,974,386 | 4,974,386 |
| **Mean autosome coverage** | 0.0773 | 0.0785 |
| **Mean Z coverage** | 0.0761 | 0.0422 | 
| **Mean Z:A ratio** | 0.98 | 0.54 | 
| **Genetic sex** | M | F |



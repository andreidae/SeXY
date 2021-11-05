# SeXY

SeXY is a simple sex-identification method for taxa lacking assembled conspecific sex chromosomes, applicaple to shotgun data from any species with a heterogametic and a homogametric sex. 
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

G. Calculate the ratio of X-linked scaffold mean coverage to Autosomal scaffold mean coverage (X:A).

H. Repeat coverage calculations (E-G) ten times. 

I. Calculate the mean and standard deviations of the X:A.


## Pipeline steps

### Required software
Satsuma synteny (http://satsuma.sourceforge.net/)

BEDtools (https://bedtools.readthedocs.io/en/latest/)

SAMtools (http://samtools.sourceforge.net/)


### Additional software
BBmap (https://jgi.doe.gov/data-and-tools/bbtools/)

### Run satsuma synteny to find which scaffolds align with sex chromosomes 
`sh SatsumaSynteny.sh $1 $2 $3 $4 $5`

```
Replace each number with the following variables
$1 - Threads
$2 - Query sequence (reference genome assembly)
$3 - Target sequence (reference sex-chromosome assembly)
$4 - Output_directory
$5 - Satsuma directory
```
 - To reduce memory and time requirements you can remove all scaffolds <10kb from the reference genomes prior to alignment. This can be done easily with bbtools for example

`reformat.sh in=file.fasta out=file_10kb.fasta minlength=10000`

### Create bed files of X-scaffolds and autosomal scaffolds
- Extract regions mapping to X and Y from the satsuma output independently

`grep chromosome_Y satsuma_summary.chained.out | awk '{print $4"\t"$5"\t"$6}' > Y.bed`

`grep chromosome_X satsuma_summary.chained.out | awk '{print $4"\t"$5"\t"$6}' > X.bed`

- Remove regions that overlapping regions from the X bed file (putatively pseudoautosomal regions)

`grep -v -f <(bedtools intersect -a X.bed -b Y.bed) X.bed > X_trim.bed`

- Extract scaffolds not mapping to X and Y from the satsuma output independently
 
`samtools faidx reference.fasta`

`grep -v -f <(cat Y.bed X.bed | cut -f 1 | sort | uniq) reference.fasta.fai | awk '{print $1"\t1\t"$2}'  > Autosomes.bed`


### Calculate depths of bam files

-  Randomly sample 1M sites 10x and calculate depth from said sites

`sh Coverage_calculation.sh Autosomes.bed X_trim.bed Output_directory Bamfile Output_prefix`

This will result in a txt file called `Output_directory/Outputprefix_ratios.txt` where the 10 X:A ratios can be found




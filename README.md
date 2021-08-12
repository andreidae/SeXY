# SeXY

SeXY is a simple sex-identification method for taxa lacking assembled conspecific sex chromosomes, applicaple to shotgun data from any species with a heterogametic and a homogametric sex. 
It uses a synteny-based approach to identify putative X-linked scaffolds in a given reference assembly, and take advantage of the expectation that males have half the amount of X-chromosome genetic material compared to females.

The seXY method requires: 
1. Raw sequencing reads of the target individual, 
2. an assembled genome of the target species/or closely related species hereafter referred to as the ‘reference genome assembly’, and 
3. assembled X and Y chromosomes from a related chromosome-level assembly hereafter referred to as the ‘reference sex-chromosome assembly’.



In general, the steps of seXY are:

A. Use the reference sex-chromosome assemblies (X and Y) to identify putative sex-linked scaffolds in the reference genome assembly via synteny.
B. Extract scaffolds aligning to the X reference sex chromosome assembly and remove pseudoautosomal regions, i.e. short homology regions between the X and Y chromosome. 
C. Extract autosomal scaffolds, i.e. scaffolds not aligning to either the X or Y reference sex-chromosome assembly. 
D. Map raw shotgun reads to the reference genome assembly.
E. Calculate the mean coverage of ten million sites randomly sampled across the extracted X-linked scaffolds. 
F. Calculate the mean coverage of ten million sites randomly sampled across the extracted autosomal scaffold.
G. Calculate the ratio of X-linked scaffold mean coverage to Autosomal scaffold mean coverage (X:A).
H. Repeat coverage calculations (E-G) ten times. 
I. Calculate the mean and standard deviations of the X:A.


# raccoon - A scalable reference assembly error correction pipeline

***THIS IS AN EARLY AND BUGGY ALPHA RELEASE NOT INTENDED FOR DISTRIBUTION!***

### What does raccoon do?

Raccoon is a pipeline that is intended to deployed as the last step of a reference genome assembly ('polishing'). It takes as input the raw data used to assemble the reference, as well as the assembly. The reads used to perform variant calls of single nucleotide variants (SNV) and short (<7bp) insertions and deletions (indels). After hard filtering, the called variants are integrated into the reference and the raw data is remaped onto the modifed assembly. Decreases in observed read depths of full length matching reads in regions were a variant was integrated are used to discard wrongly called variants. Error correction of reference assemblies is currently particularly important for regions where PacBio data has been used, as these contain an elevated amount of indels, even after self correction. 
Racoon is implemented to drive and 
This algorythm has been described by 

### Instalation

Racoon is written in python 3 (tested with version 3.4). It currently depends on the follwing external programs (parenthesis denote the version used in development. This are the versions you should use, any other ones are not guaranteed to work):

- samtools (v1.2)
- bwa (v0.7.8)
- picardtools (v.1.136)
- gatk (v3.4)

Furthermore, the follwing non-default python3 modules are required:

- pysam (v.0.8.3)
- pyvcf (v.0.6.7)

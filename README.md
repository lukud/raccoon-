# raccoon - A Scalable Reference Assembly Error Correction Pipeline

***THIS IS AN EARLY AND BUGGY ALPHA RELEASE NOT INTENDED FOR DISTRIBUTION!***


### Known bugs and To-Dos:
- Number of jobs in cluster submission may not be smaller then number of scaffolds in the assembly
- Neither number of jobs nor number scaffolds may be one
- Needs to be able to parse various lanes/readfiles in same protocol
- Need to organize all output of each stage in folders (it's rather messy right now)
- Organize module structure and separate from driver
- Need's a lot more input validation
- A deamon to automate submission of scattered stages would be nice
- do a dedup stage in first pass and extract reads from that for the follwing one?


### What does raccoon do?

Raccoon is a semi automated pipeline that is intended to be deployed as the last step of a reference genome assembly in order to correct base substitutions and small indels ('polishing'). It takes as input the raw data used to assemble the reference, as well as the assembly itself. The reads used to perform variant calls of single nucleotide variants and short (<7bp) indels. After hard filtering, the called variants are integrated into the reference and the raw data is remapped onto the modifed assembly. Decreases in observed read depths of full length matching reads in regions were a variant has been integrated are used as a criterium to reject wronlgy performed error corrections. This process is performed iteratively.
Polishing reference assemblies is currently particularly important for regions where PacBio data has been used, as these contain an elevated amount of indels even after self correction with quiver. 

### More details..

Raccoon is implemeted as a driver with accessory scripts, that executes or automatically submits to a distributed cluster. One pass of error corrections is comprised of 12 stages, with one additionall setup stage that only needs to be executed once. Wherever possible, the stages have been parallelized via a scatter-gather approach to make error corrections of large genomes feasible. Each stage automatically executes the follwing one until it hits a scattered stage. These scattered stages are highlited in yellow in the follwing schmee. The user only needs to execute stages follwing a scattered stage.

![alt tag](./pics/raccoon-scheme.png)


The basic algorithm upon which raccoon is based has first been described by [Otto et al., 2010](http://www.ncbi.nlm.nih.gov/pubmed/20562415). However, its current implementation ([iCORN2](http://icorn.sourceforge.net)) fails to scale for genomes beyond 300-400Mb, a limitation that we seek to overcome with raccoon.

### Instalation

Racoon is written in python 3 (tested with version 3.4). It currently depends on the follwing external programs (parenthesis denote the version used in development. This are the versions you should use, any other ones are not guaranteed to work):

- [samtools (v1.2)](http://www.htslib.org/download/)
- [bwa (v0.7.8)](http://sourceforge.net/projects/bio-bwa/files/)
- [picardtools (v.1.136)](http://broadinstitute.github.io/picard/)
- [gatk (v3.4)](https://www.broadinstitute.org/gatk/download/)

Furthermore, the follwing non-default python3 modules are required:

- [pysam (v.0.8.3)](https://pypi.python.org/pypi/pysam)
- [pyvcf (v.0.6.7)](https://pypi.python.org/pypi/PyVCF)

For the actual pipeline, there is currently no installation process. Just download the folder and export it to you `$PATH` like so:
```
cd /path/to/raccoon
export PATH=$PATH:$(pwd)
```

### Running raccoon

The pipeline is run through the driver like so:
```
raccoon stage protocol
```
The stage argument is one of the follwing 12:
- setup
- index
- *map*
- merge
- prepvarcall
- *varcall*
- varintegration
- reindex
- *remap*
- remerge
- correction

The setup stage only needs to be run once at the very beginning. The remaining stages can be run iteratively, raccoon will automatically take care of rewiring the input for each iteration. By default, each stage will automatically call the following stages once it finishes, until a scattered stage is reacher (denoted in bold letters above). This means that calling index will automatically call index and map. Calling map will call merge will call merge, prepvarcall and varcall. Calling varintegration will call varintegration, reindex and remap. Calling remerge will call remerge, correction and subsequently index and map for *the following iteration*. If you want to call each stage manually for some reasone (e.g. an intermediate stage failed), this can be done by invocing the -p argument, like so:
```
raccoon stage protocol -p
```


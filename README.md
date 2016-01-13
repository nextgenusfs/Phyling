PHYling
=======
Phylogentic analysis of unknown samples for the purpose of identification

Required software:
==================
* EMBOSS
* FASTX toolkit
* Muscle
* FastTree
* HMMER 3.0
* HMMER 2.0
* python
* perl - BioPerl modules
* cdbfasta
* genewise/genewisedb
* Phrap

PIPELINE FLOW:
==============
This pipeline takes an unknown sample in either fasta or fastq (_with a phred score for conversion_)
and converts the DNA sequences to 6 frame translated proteins for searching the HMMs of the markers
each organims in the databse. Contigs are created from the top hits and searched against databse
again. The top hit contig is then combined with a file of the marker protein sequences and aligned.
The alignment is then built in to a tree for viewing and comparison. There is also a script that 
will return a text document of the nearest neighbor for the unknown.

Usage:
======
For the phyling pipeline
```
./phyling     # give the usage message
./phyling <input.fasta>
./phyling <input.fastq> <phred>
```
For the marker update pipeline
```
./marker_update.sh <path to marker files>
```

Folder Dependencies
==================
Each of these scripts (in the current state) expect the following file structure:
```
Main folder
    marker_DB/
        hmm2/
        hmm3/
    results/
    databases/
        Fungal_genomes/
```

This can all be changed with a find and replace to suit your system. Alternatively symbolic links
to the appropriate files could also be created in the main folder to allow the same functionality.

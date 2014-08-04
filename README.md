PHYling
=======
Phylogentic analysis of unknown samples for the purpose of identification

Required software:
==================
*EMBOSS\n
*FASTX toolkit\n
*Muscle\n
*FastTree\n
*HMMER 3.0\n
*HMMER 2.0\n
*python\n
*perl - BioPerl modules\n
*cdbfasta\n
*genewise/genewisedb\n
*Phrap\n

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
For the phyling pipeline\n
```
./phyling -> give the usage message\n
./phyling <input.fasta>\n
./phyling <input.fastq> <phred>\n
```
For the marker update pipeline\n
``
./marker_update.sh <path to marker files>\n
``

Folder Dependencies
==================
Each of these scripts (in the current state) expect the following file structure:\n
Main folder\n
1. marker_DB/
	1. hmm2/
	2. hmm3/
2.results/
3.databases/Fungal_genomes/

This can all be changed with a find and replace to suit your system. Alternatively symoblic links
to the appropriate files could also be created in the main folder to allow the same functionality.

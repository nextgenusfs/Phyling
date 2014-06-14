#! /bin/bash
##V5.4 - updated for fastq files of unknown origin
######################################
#Stephen Bolaris 
#UC Riverside GGB - Bioinformatics
#pipeline for phylogentic analysis via 
#protein marker comparison to determine
#identitity of unknown organism
########################################

file_path="$1"
phred="$2"

#tbl="../results/hmm_matches/markers-vs-Anid.tbl"
if [ $# -eq 0 ]
then
echo "Input file path required"
echo "Usage: ./phyling.sh <input.fasta>"
echo "usage: ./phyling.sh <input.fastq> <phred>"
exit 1
fi

file_name=${file_path##*/}
base_name=${file_name%.*}

ext=${file_name##*.}
echo $ext
echo $phred
if [ "$ext" == "fastq" ]
then
#fastq file processing covert to fasta format
module load fastx_toolkit
fastq_to_fasta -Q $phred -i $file_path  > $base_name.fasta
file_path=$base_name.fasta
file_name=$base_name.fasta
fi

#This script does 3 things no names get replaced with unknonwn+#, adds UNK| to all Ids for 
#building the consesnus tree, and removes :'s that cause issues for emboss
python name_checker.py $base_name.fasta

mv $base_name.fx.fasta $base_name.fasta


#get base ID for later on when searching for hmm matches
id=`grep -o '>'[A-Z,a-z,0-9]* $file_name|sed "s/>//g"|sort -u`
echo $id
#tranlate in 6 frames the sequnences in the files
module load EMBOSS
transeq -frame 6 $file_name  $base_name.6frame.pep
# get the names from your HMM search results - only want from Anid1 for now so restricting to # those with IYBSL in the name
module load hmmer/3.0
`hmmsearch -E 1e-3  --domtblout $base_name.MARKERS.tbl ../marker_DB/markers_3.hmm $base_name.6frame.pep > $base_name.hmmsearch.out`
tbl="$base_name.MARKERS.tbl"
module load cdbfasta
# index the db for retrieval
cdbfasta $file_name
#for each marker 
while read marker
do
#get the top hits for each marker
grep $marker $tbl | grep $id | awk '{print $1}' | perl -p -e 's/_[1-6]$//'|sort -u > $base_name.$marker.read_names.lst
#get the get the sequences of the matches
cat $base_name.$marker.read_names.lst | cdbyank $file_name.cidx  > $base_name.$marker.r1.fasta
#build contigs (if any) for best results
/opt/phred-phrap-consed/bin/phrap $base_name.$marker.r1.fasta
#do a genewise search with the hmm for best peptide 
module load wise
genewisedb -hmmer -pep ../marker_DB/hmm2/$marker.hmm $base_name.$marker.r1.fasta.contigs > $base_name.$marker.candidate.pep 
#align with markerfile
cat ../marker_DB/marker_files/$marker.fa >$base_name.$marker.fasta

#get top hit from genewise peptides
awk "/>$base_name/,/\/\//" $base_name.$marker.candidate.pep >$base_name.$marker.pep
#add new line to local marker file to make sure peptide will be on it ownline
echo "\n" >> $base_name.$marker.fasta
#add top genewise peptide to local marker file
while read line; do  if [ "$line" == "//" ]; then  break; fi; echo $line >> $base_name.$marker.1.pep; done < $base_name.$marker.pep
sed "s/>/>UNK|/g" $base_name.marker.1.pep
cat $base_name.$marker.1.pep >> $base_name.$marker.fasta
#build tree
load module FastTree
load module muscle
#align local marker file with new peptide from genewise
muscle -in $base_name.$marker.fasta -out $base_name.$marker.msa
#create tree for specific marker
FastTree $base_name.$marker.msa > $base_name.$marker.tree
#get name of contig to nearest neighbor search 
contig=`grep '>' $basen_name.$marker.1.pep| sed "s/>//g"`
neigh=`perl get_nearest_tree_neighbor.pl $base_name.$marker.tree $contig`
echo "$base_name.$marker.tree\t$neigh\n" >> $base_name.closest_match.txt

done < marker_list.txt
#done < marker_test.txt

#clean usage directory
stamp=`date +%D|sed "s/\//_/g"`
mkdir $base_name\_$stamp
mv $base_name*\.* $base_name\_$stamp
mv $base_name\_$stamp ../results/


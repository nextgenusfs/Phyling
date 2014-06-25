#/bin/bash
#Stephen Bolaris v1.0

CUTOFF=1e-100
DB="../DB/genomes/fungi"
LIST=fungi_genome_list.txt
marker_path=$1
#echo $marker_path
#read -p "check, check"

if [ $# -ne 1 ]
then
 echo "Path to markers required as argument to script\n"
 exit 1
fi

#
# cp marker files in to pwd -> forward 
# could just append to that file by passing the path to the python script 
# as a second paramater
while read line
do
# cp $marker_path/$line/*.fasta ./$line.fasta
 cat $marker_path/$line/*.fasta > $line.fasta
done < marker_list.txt

#for each of the 164 files
#ls $DB > genome_list.txt

#load required modules
module load hmmer/3.0
module load hmmer/2.3.2
module load cdbfasta
module load muscle

for file in $DB
do

base_name=${genome%.*}
#get that file and move to pwd
cp $db_loc$genome /tmp
#make the fasta file in to a DB to yank sequence
cdbfasta /tmp/$genome
#hmm search that file to find the peptide with highest match
# for each marker, using the markers.hmm
hmmsearch -E $CUTOFF --domtblout $base_name.MARKERS.tbl ../marker_DB/markers_3.hmm /tmp/$genome > $base_name.hmmsearch.out

#get the peptide that matches each marker
#make a python script to do this
python hmm_marker_update.py $base_name.MARKERS.tbl
#add (append) the new peptide to the new marker file
#clean up the genome file (ie remove the copy that is script folder)
rm $base_name.*
done
# < genome_list.txt

#provides a double check that duplicate seqeunces are not added
while read marker
do
perl uniq_marker.pl $marker.fasta
done < marker_list.txt

#move currrent maker files to a backup location with date
#move new markerfiles to thier location
DATE=`date +%d_%m_%y`
mv ../marker_DB/hmm2 ../marker_DB/hmm2_$DATE
mv ../marker_DB/hmm3 ../marker_DB/hmm3_$DATE 
#make directories for where the two different types of HMMs will go
mkdir ../marker_DB/hmm2
mkdir ../marker_DB/hmm3
#depricated for debugging purposes
#read -p "check marker files"

#for each marker in the list
while read marker
do
 #create MSA with muscle
 muscle -in $marker.fasta -out $marker.msa

 #create new HMMs (v3 and v2)
 hmmbuild2 $marker.hmm $marker.msa

 mv *.hmm ../marker_DB/hmm2/
 hmmbuild --informat afa  $marker.hmm $marker.msa
 mv *.hmm ../marker_DB/hmm3/

 #move old HMMs to back up location with date
 #move new HMMs in to thier spot
done < marker_list.txt

#read -p "check movement"
mkdir ../marker_DB/marker_update_$DATE
mv *.msa ../marker_DB/marker_update_$DATE

#clean up workspace
while read line
do
 rm /tmp/$line*
done < genome_list.txt

rm ../marker_DB/*.hmm
cat ../marker_DB/hmm2/*.hmm > ../marker_DB/markers_2.hmm
cat ../marker_DB/hmm3/*.hmm > ../marker_DB/markers_3.hmm

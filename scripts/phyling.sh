#! /bin/bash
##V5.4 - updated for fastq files of unknown origin
######################################
#Stephen Bolaris 
#UC Riverside GGB - Bioinformatics
#pipeline for phylogentic analysis via 
#protein marker comparison to determine
#identitity of unknown organism
########################################
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DBDIR=`dirname $BASEDIR/`"/DB"
CLADE=fungi
CUTOFF=1e-3
PHRAP=/opt/phred-phrap-consed/bin/phrap
# modules loading for phyling
module load EMBOSS
module load hmmer/3.0
module load cdbfasta
module load fastx_toolkit
module load wise
module load FastTree
module load muscle
module load trimal

file_path="$1"
phred="$2"


if [ $# -eq 0 ]; then
 echo "Input file path required"
 echo "Usage: ./phyling.sh <input.fasta>"
 echo "usage: ./phyling.sh <input.fastq> <phred>"
 exit 1
fi

file_name=${file_path##*/}
base_name=${file_name%.*}

ext=${file_name##*.}
#echo "extension is $ext"
#echo "phred is $phred"

if [ "$ext" == "fastq" ]; then
 #fastq file processing covert to fasta format
 if [ ! -f $base_name.fa ]; then
  fastq_to_fasta -Q $phred -i $file_path  > $base_name.fasta
 fi
 file_name=$base_name.fasta
fi

#This script does 3 things no names get replaced with unknown+#
# adds UNK| to all Ids for building the consesnus tree
# and removes :'s that cause issues for emboss
if [ ! -f $base_name.fx.ok ]; then
 perl $BASEDIR/name_checker.pl $file_name > $base_name.fx
 mv $base_name.fx $file_name
 touch $base_name.fx.ok
fi

#get base ID for later on when searching for hmm matches
id=`grep -o '>'[A-Z,a-z,0-9]* $file_name|sed "s/>//g"|sort -u`
#echo $id

#tranlate in 6 frames the sequnences in the files
if [ ! -f $base_name.6frame.faa ]; then
 transeq -trim -clean -frame 6 $file_name $base_name.6frame.faa
 perl -i -p -e "s/^>/>$id\|/" $base_name.6frame.faa
fi
# get the names from your HMM search results - only want from Anid1 for now so restricting to # those with IYBSL in the name

tbl="$base_name.MARKERS.tbl"
if [ ! -f $tbl ] || [ ! -s $tbl ]; then
 hmmsearch -E $CUTOFF  --domtblout $tbl \
  $DBDIR/markers/$CLADE/markers_3.hmmb $base_name.6frame.faa \
  > $base_name.hmmsearch.out
fi


# index the db for retrieval
if [ ! -f $file_name.cidx ]; then
 cdbfasta $file_name
fi

#for each marker
while read marker; do
 #get the top hits for each marker
 grep $marker $tbl | grep $id | awk '{print $1}' | perl -p -e 's/_[1-6]$//'|sort -u > $base_name.$marker.read_names.lst
#get the get the sequences of the matches
cat $base_name.$marker.read_names.lst | cdbyank $file_name.cidx  > $base_name.$marker.r1.fasta

 if [ -s $base_name.$marker.r1.fasta ]; then
  #build contigs (if any) for best results
  if [ ! -f $base_name.$marker.r1.fasta.contigs ]; then
   $PHRAP $base_name.$marker.r1.fasta
  fi
 fi

 if [ -s $base_name.$marker.r1.fasta.contigs ]; then
  #do a genewise search with the hmm for best peptide 
   if [ ! -f  $base_name.$marker.candidate.pep ]; then
    genewisedb -hmmer -pep -splice flat -init local -silent -quiet $DBDIR/markers/$CLADE/HMM2/$marker.hmm $base_name.$marker.r1.fasta.contigs > $base_name.$marker.candidate.pep 
   fi
 fi

 if [ -s $base_name.$marker.candidate.pep ]; then 
  #get top hit from genewise peptides
  awk "/>$base_name/,/\/\//" $base_name.$marker.candidate.pep > $base_name.$marker.pep

 # this could be re-written as a cleaner script (e.g. perl/python: genewise2pep.pl)
 #add top genewise peptide to local marker file
 while read line; do  
  if [ "$line" == "//" ]; then 
   break; 
  fi
  echo $line >> $base_name.$marker.1.pep
 done < $base_name.$marker.pep

 fi

 if [ -s $base_name.$marker.1.pep ]; then
  #align with markerfile
  cp $DBDIR/markers/$CLADE/marker_files/$marker.fa $base_name.$marker.fasta
  #add new line to local marker file to make sure peptide will be on it ownline
  echo "" >> $base_name.$marker.fasta
  
  perl -i -p -e "s/>/>$id|/g" $base_name.$marker.1.pep 
  cat $base_name.$marker.1.pep >> $base_name.$marker.fasta

  #build tree
  # MSA
  #align local marker file with new peptide from genewise
  # muscle -in $base_name.$marker.fasta -out $base_name.$marker.msa
  if [ ! -f $base_name.$marker.msa ]; then 
   hmmalign --trim --amino $DBDIR/markers/$CLADE/HMM3/$marker.hmm $base_name.$marker.fasta > $base_name.$marker.msa
   sreformat clustal $base_name.$marker.msa > $base_name.$marker.aln
   mv $base_name.$marker.aln $base_name.$marker.msa
  fi
  if [ ! -f $base_name.$marker.msa.trim ]; then
   trimal -in $base_name.$marker.msa -out $base_name.$marker.msa.trim -automated1 -fasta
  fi
 #create tree for specific marker
 # FastTree $base_name.$marker.msa > $base_name.$marker.tree
  #FastTree -boot 1000 -wag -bionj $base_name.$marker.msa.trim > $base_name.$marker.tree
 #get name of contig to nearest neighbor search 
 # contig=`grep '>' $basen_name.$marker.1.pep | sed "s/>//g"`
 # neigh=`perl $BASEDIR/get_nearest_tree_neighbor.pl $base_name.$marker.tree $contig`
 # perl -e "print \"$base_name.$marker.tree\t$neigh\n\";" >> $base_name.closest_match.txt
 fi # deal with no results more swiftly

#DEBUG
# if [ -f $base_name.$marker.fasta ]; then
#  break
# fi

done < $DBDIR/markers/$CLADE/list.txt

perl $BASEDIR/combine_fasaln.pl -d . -ext msa.trim -if fasta -of fasta -o combine.all.fasaln

#clean usage directory
#stamp=`date +%D|sed "s/\//_/g"`
#mkdir $base_name\_$stamp
#mv $base_name*\.* $base_name\_$stamp
#mv $base_name\_$stamp ../results/


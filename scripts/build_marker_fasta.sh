#PBS -l nodes=1:ppn=4 -N hmm_rebuild_markers -j oe
CPU=1
F=$PBS_ARRAYID
CUTOFF=1e-10
MARKERS=markers_3.hmmb
MARKERDIR=DB/markers/fungi
GENOMEDIR=DB/genomes/fungi
GENOMELIST=$GENOMEDIR/fungi_genome_list.txt
OUT=out
if [ ! -d $OUT ]; then
 mkdir -p $OUT
fi

if [ ! $F ]; then
 F=$1
fi

if [ ! $F ]; then
 echo "no PBS_ARRAYID or input"
 exit
fi

if [ $PBS_NUM_PPN ]; then
 CPU=$PBS_NUM_PPN
fi

FILE=`head -n $F $GENOMELIST | tail -n 1`
if [ ! $FILE ]; then 
 echo "No input file - check PBS_ARRAYID or input number"
 exit
fi

BASE=`basename $FILE .fasta`
echo "FILE is $FILE. Base is $BASE"

if [ ! -f "$OUT/$BASE.domtbl" ]; then
 hmmsearch --domtbl $OUT/$BASE.domtbl -E $CUTOFF --cpu $CPU $MARKERDIR/$MARKERS $GENOMEDIR/$FILE > $OUT/$BASE.hmmsearch
fi

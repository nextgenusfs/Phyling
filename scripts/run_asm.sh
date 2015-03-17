#PBS -l nodes=1:ppn=2 -N phyling -j oe 
module load perl
FOLDER=input_asm
N=$PBS_ARRAYID
INFILE=lst_asm
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "Need a PBS_ARRAYID or cmdline number"
 exit;
fi

line=`head -n $N $INFILE  | tail -n 1`

perl ../../Phyling/scripts/phyling.pl $FOLDER/$line --cpu $PBS_NP

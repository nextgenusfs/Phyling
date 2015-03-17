#PBS -l nodes=1:ppn=4 -N tcoffee -j oe
F=$PBS_ARRAYID
GENOMELIST=FILES
if [ ! -f $GENOMELIST ]; then
 ls *.aa.rename > $GENOMELIST
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

t_coffee $FILE -n_core $CPU

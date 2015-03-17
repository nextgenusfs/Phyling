#PBS -l nodes=1:ppn=24 -q js -N raxmlPep -j oe
module load RAxML/8.1.1

CPU=2

if [ $PBS_NUM_PPN ]; then
 CPU=$PBS_NUM_PPN
fi

raxmlHPC-PTHREADS-SSE3 -T $CPU -# 100 -x 121 -f a -p 123 -m PROTGAMMAAUTO -s allseq_pep.nogap.phy -n allseq_pep_nogap_ML_GAMMA -o Gonapodya_prolifera,Catenaria_anguillulae
raxmlHPC-PTHREADS-SSE3 -T $CPU -# 100 -x 121 -f a -p 123 -m PROTGAMMAAUTO -s allseq_pep.strictplus.phy -n allseq_pep_strictplus_ML_GAMMA -o Gonapodya_prolifera,Catenaria_anguillulae

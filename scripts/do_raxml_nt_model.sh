#PBS -l nodes=1:ppn=24 -q js -N raxmlCodonModel -j oe
module load RAxML/8.1.1

CPU=2

if [ $PBS_NUM_PPN ]; then
 CPU=$PBS_NUM_PPN
fi

raxmlHPC-PTHREADS-SSE3 -T $CPU -# 100 -x 121 -f a -p 123 -m GTRCATX -s allseq.strictplus.phy -n allseq_strictplus_codonmodel_ML_Codons -q codon_model_strictplus.txt -B 0.03 -o Gonapodya_prolifera,Catenaria_anguillulae
raxmlHPC-PTHREADS-SSE3 -T $CPU -# 100 -x 121 -f a -p 123 -m GTRCATX -s allseq.nogap.phy -n allseq_strictplus_codonmodel_ML_Codons -q codon_model_strictplus.txt -B 0.03 -o Gonapodya_prolifera,Catenaria_anguillulae

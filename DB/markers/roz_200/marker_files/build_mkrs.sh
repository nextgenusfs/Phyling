#PBS -j oe -N buildMrkPhylingRoz

perl ../../../../scripts/aln2hmm.pl -g roz_200 -i .

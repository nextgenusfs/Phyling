#PBS -j oe -N buildMrkPhylingRoz

perl ../../../../scripts/aln2hmm.pl -g AFTOL1 -i .

#PBS -j oe -l walltime=2:00:00
perl ../../Phyling/scripts/gather_besthmm_build_markerfasta.pl -pep outgroups/pep -o outgroup_markers -i outgroups -cds outgroups/cds

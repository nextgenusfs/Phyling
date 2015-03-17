
cd marker_aln
grep -h ">" *.aa.rename | awk '{print $1}' | perl -p -e 's/^>//' | sort | uniq > ../strains

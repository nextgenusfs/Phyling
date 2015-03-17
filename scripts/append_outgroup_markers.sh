for file in outgroup_markers/*.aa.rename; do 
 x=`basename $file .aa.rename`
 cat $file >> marker_aln/$x.aa.rename
done

for file in outgroup_markers/*.cds.rename; do
 x=`basename $file .cds.rename`
 cat $file >> marker_aln/$x.cds.rename
done


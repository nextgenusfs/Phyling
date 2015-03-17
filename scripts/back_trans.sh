module load trimal
for file in *.aa.aln; do 
x=`basename $file .aa.aln`;
echo $file
trimal -in $file -out $x.cds.fasaln_strictplus.trim -splitbystopcodon -backtrans $x.cds.rename -strictplus -fasta
trimal -in $file -out $x.cds.fasaln_strict.trim -splitbystopcodon -backtrans $x.cds.rename -strict -fasta
trimal -in $file -out $x.cds.fasaln_nogap.trim -splitbystopcodon -backtrans $x.cds.rename -nogaps -fasta
done

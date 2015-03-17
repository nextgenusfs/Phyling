module load trimal
for file in *.aa.aln; do 
x=`basename $file .aa.aln`;
trimal -in $file -out $x.fasaln_strictplus.pep.trim -strictplus -fasta
trimal -in $file -out $x.fasaln_strict.pep.trim -strict -fasta
trimal -in $file -out $x.fasaln_nogap.pep.trim -nogaps -fasta
done

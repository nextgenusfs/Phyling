module load trimal
for file in *.aln
do 
 trimal -in $file -out $file.trim -fasta -automated1
done

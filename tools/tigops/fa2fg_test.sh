./tigops fasta2fastg  -tigs test/graph.fa  -kmer-size 3 -out test/out.fastg
echo "input"
cat test/graph.fa
echo "output"
cat test/out.fastg

./tigops fasta2fastg  -tigs test/X.fa -kmer-size 4 -out test/out.fastg
echo "input"
cat test/X.fa
echo "output"
cat test/out.fastg

input=test/Xrev.fa
./tigops fasta2fastg  -tigs $input -kmer-size 4 -out test/out.fastg
#echo "input"
#cat $input
echo "output"
cat test/out.fastg

./tigops fasta2fastg  -tigs test/ecoli_spades_contigs_k33.fasta  -kmer-size 34 -out test/out.fastg # setting k=34 because spades use edge-centric dbg
grep ">" test/out.fastg > test/out.fastg.ids
diff test/out.fastg.ids test/ecoli_spades_contigs_k33.fastg.ids

var=$?

if [ $var -eq 0 ] 
then
    echo Test PASSED
    exit 0
else
    echo Test FAILED
    exit 1
fi

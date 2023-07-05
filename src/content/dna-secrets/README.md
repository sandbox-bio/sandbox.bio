## Notes

```bash
# Encode string
perl dna-encode.pl <(echo "sandbox")

# TC  TCTCATTCTG CTGCC TGCCACCGACGGTTCGTCCATTTCGGC
# CT  CTCTTAGTTC CTGCC GCCCTGCCATCCGTATGTTCTGGGCCA

# CTGGCATTATTAACCTCTTAGTTCCTGCCGCCCTGCCATCCGTATGTTCTGGGCCA
# C G C T A T A C C T G T C G C C C G C T C T T T C G G C
#  T G A T T A C T T A T C T C G C T C A C G A G T T G C A

# Get lambda ref
cp ~/Documents/dev/wasm/biowasm/tools/bowtie2/src/example/reference/lambda_virus.fa ref1.fa
cp ~/Documents/dev/wasm/biowasm/tools/bowtie2/src/example/reference/lambda_virus.fa ref2.fa
# split refs into 35-chars lines

# Modify ref to introduce SNPs
vi ref1.fa  # CGCTATACCTGTCGCCCGCTCTTTCGGC  # even line numbers
vi ref2.fa  # TGATTACTTATCTCGCTCACGAGTTGCA  # odd line numbers

# Simulate reads; easiest = no mutations :)
wgsim -N 300 -e 0 -r 0 -R 0 -X 0 -S 42 -1 40 ref1.fa reads1.fq /dev/null
wgsim -N 300 -e 0 -r 0 -R 0 -X 0 -S 42 -1 40 ref2.fa reads2.fq /dev/null
cp reads1.fq ~/Documents/dev/sandbox.bio/src/tutorials/dna-secrets/data/reads.fq
cp reads2.fq ~/Documents/dev/sandbox.bio/src/tutorials/dna-secrets/data/morereads.fq
./build.sh
```

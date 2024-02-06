## Notes

```bash
# Download data
curl -o data.fasta "https://raw.githubusercontent.com/osuoer/computational-biology/master/orf_trans.fasta"

# Subsample
head -n 500 data.fasta > orf_trans.fasta

# Add sequences needed for exercise
seqtk subseq -l 60 data.fasta <(echo -e "YHR007C\nYDR402C") >> orf_trans.fasta
```

## Notes

```bash
# Subset
samtools view -h original.bam 20:1.3e6-1.5e6 > sample.sam

# Shuffle it so it's realistic
cat <(samtools view -H sample.sam) <(shuf <(samtools view sample.sam)) > sample.shuffled.sam
mv sample.shuffled.sam sample.sam

# But want to keep some of the non-properly-paired reads
samtools view -F 0x2 original.bam | head -n5 >> sample.sam
samtools view -b sample.sam | samtools sort > sample.bam
samtools index sample.bam

## Downsample (but results in very low coverage...)
##samtools view -h -S -s 42.0001 ~/Downloads/sample.sam > sample.small.sam
```

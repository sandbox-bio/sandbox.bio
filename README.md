# sandbox.bio

## Tutorials

### bedtools-tutorial

http://quinlanlab.org/tutorials/bedtools/bedtools.html

```bash
# Download data
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/maurano.dnaseI.tgz
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/cpg.bed
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/exons.bed
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/gwas.bed
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/genome.txt
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/hesc.chromHmm.bed

# Subsample
for file in $(ls -f .); do echo $file; head -n200 $file > downsampled/$file; done
echo -e "chr1\t30365652" > downsampled/genome.txt

cp cpg.bed exons.bed gwas.bed hesc.chromHmm.bed genome.txt fHeart-DS16621.hotspot.twopass.fdr0.05.merge.bed fHeart-DS15839.hotspot.twopass.fdr0.05.merge.bed fSkin_fibro_bicep_R-DS19745.hg19.hotspot.twopass.fdr0.05.merge.bed ../../tutorials/1-intro-to-bedtools/data
```

#### What's differerent
* Removed Setup (files are preloaded for users)
* Using 3/20 of the bed files + Downsampled them to 200 lines (and genome.txt to just chr1) to make it easier for me :)

---

### bowtie2-tutorial

```bash
REF=/bowtie2/example/index/lambda_virus
bowtie2 -x $REF -U reads_1.fq -S eg1.sam
```

#### What's differerent
* Removed Setup (files are preloaded for users)
* Subsampled reads_1.fq, reads_2.fq and longreads.fq

---

### samtools tutorial

```bash
# Subset
samtools view -h original.bam 20:1.3e6-1.5e6 > sample.sam
samtools view -b -o sample.bam sample.sam
samtools index sample.bam

# Shuffle it so it's realistic
cat <(samtools view -H sample.sam) <(shuf <(samtools view sample.sam)) > sample.shuffled.sam

# Downsample (but results in very low coverage...)
samtools view -h -S -s 42.0001 ~/Downloads/sample.sam > sample.small.sam
```

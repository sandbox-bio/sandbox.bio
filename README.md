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

TODO: Add `bcftools` to tutorial:

```bash
bowtie2  -x $REF  -1 reads_1.fq  -2 reads_2.fq  -S eg2.sam; samtools view eg2.sam -o eg2.bam; samtools sort eg2.sam -o eg2.sorted.bam

bcftools mpileup -f $REF_FA eg2.sorted.bam
```

Setup:

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

---

### hidden-message

**Prepare**:

```bash
# Get lambda ref
cp ~/Documents/dev/wasm/biowasm/tools/bowtie2/src/example/reference/lambda_virus.fa ref.fa

# Modify ref to introduce SNPs
vi ref.fa

# Simulate reads; easiest = no mutations :)
wgsim -N 500 -e 0 -r 0 -R 0 -X 0 -S 42 ref.fa r1.fq r2.fq
cp r1.fq ~/Documents/dev/sandbox.bio/src/tutorials/hidden-message/data/reads.fq
cd /Users/robert/Documents/dev/sandbox.bio/ && ./build.sh && cd -

# 
```

**Solution**:

```bash
bowtie2 -x $REF -U reads.fq -S aligned.sam
samtools sort -o aligned.sorted.bam aligned.sam
samtools index aligned.sorted.bam
bcftools mpileup -f $REF_FASTA -o variants.vcf eg1.sorted.bam
bcftools call -mv -Ob -o variants.bcf variants.vcf
bcftools query -f'%ALT' variants.bcf -o secret  # output to screen doesn't work b/c not flushed
cat secret

# One liner:
bowtie2  -x $REF  -U reads.fq  -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam; samtools index aligned.sorted.bam; bcftools mpileup -f $REF_FASTA -o variants.vcf aligned.sorted.bam; bcftools call -mv -Ob -o variants.bcf variants.vcf; bcftools query -f'%ALT' variants.bcf -o secret; cat secret
```

**Converter**

```js
// https://stackoverflow.com/a/53247859
function binaryToString(input) {
  let bytesLeft = input;
  let result = '';

  // Check if we have some bytes left
  while (bytesLeft.length) {
    // Get the first digits
    const byte = bytesLeft.substr(0, 8);
    bytesLeft = bytesLeft.substr(8);
    result += String.fromCharCode(parseInt(byte, 2));
  }
  return result;
}

// https://science.sciencemag.org/content/337/6102/1628
binaryToString("CTCCATCACGCATGTACGACCAAT".split("").map(b => {
  if(b == "A" || b == "C") return "0";
  else return "1";
}).join(""))
```

# sandbox.bio

[![Tests](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml/badge.svg)](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml)

## Development

* Branches
  * `dev`: Development branch, pushing there deploys a preview site
  * `main`: Production branch; merge dev into main to deploy to sandbox.bio
* Add a tag when deploy to prd

## Testing

* `npm run test` will launch all the tests in headless way
* `cypress open` opens Cypress so you can see the tests inside Chrome
* `npm run test -- --spec "tests/test_tutorials.js"` runs a single file

## Debugging

### Local biowasm builds

```javascript
// Terminal.svelte
const TOOLS_DEFAULT = [{
	tool: "samtools",
	version: "1.10",
	urlPrefix: "http://localhost:12346/biowasm/tools/samtools/build/"
}];
```

### Local aioli builds

```javascript
// cli.js
	_aioli = await new Aioli(config.tools, {
		env: window.location.hostname == "localhost" ? "stg" : "prd",
		debug: window.location.hostname == "localhost",
		urlAioli: "http://localhost:12346/aioli/dist/aioli.worker.js"
	});
```

### Run command on load

```javascript
// Terminal.svelte
setTimeout(async () => {
  $xtermAddons.echo.setInput("samtools view -b /samtools/examples/toy.sam > bad.bam; samtools quickcheck bad.bam");
	$xtermAddons.echo.handleData("\r");
}, 1000)
```


## Tutorials

### Add a new one

* Create `src/tutorials/my-tutorial/config.js` + import it in `src/config.js`
* Add tests in `tests/test_tutorials.js`


## Tutorial notes

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

### dna-secrets

**Prepare**:

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

**Solution**:

```bash
bowtie2 -x $REF -U reads.fq -S aligned.sam
samtools sort -o aligned.sorted.bam aligned.sam
samtools index aligned.sorted.bam
bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools call -mv - > variants.vcf
bcftools query -f'%ALT' variants.vcf -o secret  # output to screen doesnt work because not flushed
cat secret

# One liner:
bowtie2  -x $REF  -U reads.fq  -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam; samtools index aligned.sorted.bam; bcftools mpileup -f $REF_FASTA -o variants.vcf aligned.sorted.bam; bcftools call -mv -Ob -o variants.bcf variants.vcf; bcftools query -f'%ALT' variants.bcf -o secret; cat secret
```

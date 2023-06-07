## Notes

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

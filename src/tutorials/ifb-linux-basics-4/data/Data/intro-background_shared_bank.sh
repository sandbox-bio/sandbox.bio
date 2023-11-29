#!/bin/bash

## Banks folders
mkdir -p /shared/data/bank/{bos_taurus,homo_sapiens,nr}

cd /shared/data/bank/bos_taurus
mkdir -p UMD3.1/{star-2.7.2b,fasta,bowtie2}
touch \
  UMD3.1/star-2.7.2b/{SAindex,chrLength.txt,chrName.txt,chrStart.txt,genomeParameters.txt} \
  UMD3.1/fasta/{Bos_taurus.UMD3.1.dna.toplevel_F.fa.fai,Bos_taurus.UMD3.1.dna.toplevel_F.fa} \
  UMD3.1/bowtie2/{Bos_taurus.UMD3.1.dna.toplevel_F.rev.1.bt2,Bos_taurus.UMD3.1.dna.toplevel_F.2.bt2,Bos_taurus.UMD3.1.dna.toplevel_F.rev.2.bt2,Bos_taurus.UMD3.1.dna.toplevel_F.1.bt2}

cd /shared/data/bank/homo_sapiens
mkdir -p \
  hg19/{hisat2,star-2.7.5a,fasta,bowtie2} \
  hg38/{hisat2,star-2.7.5a,fasta,bowtie2}
touch \
  hg19/hisat2/{hg19.1.ht2,hg19.2.ht2,hg19.3.ht2,hg19.4.ht2} \
  hg19/star-2.7.5a/{SAindex,chrLength.txt,chrName.txt,chrStart.txt,genomeParameters.txt} \
  hg19/fasta/{hg19.fa.fai,hg19.fa} \
  hg19/bowtie2/{hg19.1.bt2,hg19.2.bt2,hg19.rev.1.bt2,hg19.rev.2.bt2} \
  hg38/hisat2/{genome.4.ht2,genome.2.ht2,genome.3.ht2,genome.1.ht2} \
  hg38/star-2.7.5a/{SAindex,chrLength.txt,chrName.txt,chrStart.txt,genomeParameters.txt} \
  hg38/fasta/{hg38.fa,hg38.fa.fai} \
  hg38/bowtie2/{hg38.2.bt2,hg38.rev.2.bt2,hg38.rev.1.bt2,hg38.1.bt2}

cd /shared/data/bank/nr
mkdir -p nr_2018-09-28/{blast,fasta,diamond}
touch \
  nr_2018-09-28/fasta/nr.fsa \
  nr_2018-09-28/blast/{nr.01.psd,nr.01.ppi,nr.01.phd,nr.02.psd,nr.02.ppi,nr.02.phd} \
  nr_2018-09-28/diamond/{nr.dmnd,viral.protein_refseq_98.dmnd}

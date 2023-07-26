// Steps
import Step0 from "./steps/step00.md";
import Step1 from "./steps/step01.md";
import Step2 from "./steps/step02.md";
import Step3 from "./steps/step03.md";
import Step4 from "./steps/step04.md";
import Step5 from "./steps/step05.md";
import Step6 from "./steps/step06.md";
import Step7 from "./steps/step07.md";
import Step8 from "./steps/step08.md";
import Step9 from "./steps/step09.md";


export const config = {
	id: "ifb-linux-basics-3",
	pwd: "ifb-linux-basics-3",
	listed: false,
	name: "Linux Basics 3",
	subtitle: `by <a href="https://www.france-bioinformatique.fr/en/home/" target="_blank">French Institute of Bioinformatics</a>`,
	description: "Manipulating files and directories",
	tags: ["unix", "shell", "terminal"],
	tools: ["ls", "date"],
	difficulty: ["belt green", "belt blue"],
	steps: [
		{ name: "Introduction", component: Step0 },
		{ name: "Manipulating data", component: Step1 },
		{ name: "Digging through a large file", component: Step2 },
		{ name: "The begining of a file", component: Step3 },
		{ name: "The end of a file", component: Step4 },
		{ name: "Counting words and lines in a file", component: Step5 },
		{ name: "Searching patterns", component: Step6 },
		{ name: "Extracting colums", component: Step7 },
		{ name: "Some useful commands and tips", component: Step8 },
		{ name: "Congratulations", component: Step9 },
	],
	// For now, don't use `cd` in init, need to debug why it breaks the initialization script
	init: `
		mkdir -p bos_taurus homo_sapiens nr;
		PREFIX=bos_taurus/UMD3.1;
			mkdir -p $PREFIX/star-2.7.2b $PREFIX/fasta $PREFIX/bowtie2;
			touch $PREFIX/star-2.7.2b/SAindex $PREFIX/star-2.7.2b/chrLength.txt $PREFIX/star-2.7.2b/chrName.txt $PREFIX/star-2.7.2b/chrStart.txt $PREFIX/star-2.7.2b/genomeParameters.txt $PREFIX/fasta/Bos_taurus.UMD3.1.dna.toplevel_F.fa.fai $PREFIX/fasta/Bos_taurus.UMD3.1.dna.toplevel_F.fa $PREFIX/bowtie2/Bos_taurus.UMD3.1.dna.toplevel_F.rev.1.bt2 $PREFIX/bowtie2/Bos_taurus.UMD3.1.dna.toplevel_F.2.bt2 $PREFIX/bowtie2/Bos_taurus.UMD3.1.dna.toplevel_F.rev.2.bt2 $PREFIX/bowtie2/Bos_taurus.UMD3.1.dna.toplevel_F.1.bt2;

		PREFIX=homo_sapiens;
			mkdir -p $PREFIX/hg19 $PREFIX/hg19/hisat2 $PREFIX/hg19/star-2.7.5a $PREFIX/hg19/fasta $PREFIX/hg19/bowtie2 $PREFIX/hg38 $PREFIX/hg38/hisat2 $PREFIX/hg38/star-2.7.5a $PREFIX/hg38/fasta $PREFIX/hg38/bowtie2;
			touch $PREFIX/hg19/hisat2/hg19.2.ht2 $PREFIX/hg19/hisat2/hg19.1.ht2 $PREFIX/hg19/hisat2/hg19.3.ht2 $PREFIX/hg19/hisat2/hg19.4.ht2 $PREFIX/hg19/star-2.7.5a/SAindex $PREFIX/hg19/star-2.7.5a/chrLength.txt $PREFIX/hg19/star-2.7.5a/chrName.txt $PREFIX/hg19/star-2.7.5a/chrStart.txt $PREFIX/hg19/star-2.7.5a/genomeParameters.txt $PREFIX/hg19/fasta/hg19.fa.fai $PREFIX/hg19/fasta/hg19.fa $PREFIX/hg19/bowtie2/hg19.1.bt2 $PREFIX/hg19/bowtie2/hg19.2.bt2 $PREFIX/hg19/bowtie2/hg19.rev.1.bt2 $PREFIX/hg19/bowtie2/hg19.rev.2.bt2 $PREFIX/hg38/hisat2/genome.4.ht2 $PREFIX/hg38/hisat2/genome.2.ht2 $PREFIX/hg38/hisat2/genome.3.ht2 $PREFIX/hg38/hisat2/genome.1.ht2 $PREFIX/hg38/star-2.7.5a/SAindex $PREFIX/hg38/star-2.7.5a/chrLength.txt $PREFIX/hg38/star-2.7.5a/chrName.txt $PREFIX/hg38/star-2.7.5a/chrStart.txt $PREFIX/hg38/star-2.7.5a/genomeParameters.txt $PREFIX/hg38/fasta/hg38.fa $PREFIX/hg38/fasta/hg38.fa.fai $PREFIX/hg38/bowtie2/hg38.2.bt2 $PREFIX/hg38/bowtie2/hg38.rev.2.bt2 $PREFIX/hg38/bowtie2/hg38.rev.1.bt2 $PREFIX/hg38/bowtie2/hg38.1.bt2;

		PREFIX=nr/nr_2018-09-28/;
			mkdir -p $PREFIX/blast $PREFIX/fasta $PREFIX/diamond;
			touch $PREFIX/fasta/nr.fsa $PREFIX/blast/nr.01.psd $PREFIX/blast/nr.01.ppi $PREFIX/blast/nr.01.psd $PREFIX/blast/nr.01.phd $PREFIX/blast/nr.02.psd $PREFIX/blast/nr.02.ppi $PREFIX/blast/nr.02.psd $PREFIX/blast/nr.02.phd $PREFIX/diamond/nr.dmnd $PREFIX/diamond/viral.protein_refseq_98.dmnd;
		
		mkdir -p $HOME/test;
			echo "first test" > $HOME/test/first_file.txt;
			echo "second test" > $HOME/test/second_file.txt;
	`,
	files: [
		"data/ifb-linux-basics-3/MACS2.csv",
		"data/ifb-linux-basics-3/NC_009089.bed",
		"data/ifb-linux-basics-3/NC_009089.fasta",
		"data/ifb-linux-basics-3/SAOUHSC.fasta",
		"data/ifb-linux-basics-3/SAOUHSC.bed",
	],
};

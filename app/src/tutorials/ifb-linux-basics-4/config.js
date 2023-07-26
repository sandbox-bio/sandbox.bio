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


export const config = {
	id: "ifb-linux-basics-4",
	pwd: "ifb-linux-basics-4",
	listed: false,
	name: "Combine some bash commands",
	subtitle: `by <a href="https://www.france-bioinformatique.fr/en/home/" target="_blank">French Institute of Bioinformatics</a>`,
	description: "IFB Scenario 4",
	tags: ["unix", "shell", "terminal"],
	tools: ["ls", "date"],
	difficulty: ["belt blue"],
	steps: [
		{ name: "Introduction", component: Step0 },
		{ name: "Output of a Unix command", component: Step1 },
		{ name: "Input of a Unix command", component: Step2 },
		{ name: "Error stream and its redirection", component: Step3 },
		{ name: "Example input and output", component: Step4 },
		{ name: "Associate with pipe", component: Step5 },
		{ name: "More practice?", component: Step6 },
		{ name: "Summary", component: Step7 },
		{ name: "Congratulations", component: Step8 },
	],
	// For now, don't use `cd` in init, need to debug why it breaks the initialization script
	/*
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
	`,*/
	files: [
		"data/ifb-linux-basics-4/Data/NC_009089.bed",
		"data/ifb-linux-basics-4/Data/NC_009089.fasta",
		"data/ifb-linux-basics-4/Data/SAOUHSC.fasta",
		"data/ifb-linux-basics-4/Data/SAOUHSC.bed",
	],
};

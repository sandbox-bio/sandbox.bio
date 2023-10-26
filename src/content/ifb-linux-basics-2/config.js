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
	id: "ifb-linux-basics-2",
	pwd: "ifb-linux-basics-2",
	listed: false,
	name: "Manipulating files and directories",
	subtitle: `by <a href="https://www.france-bioinformatique.fr/en/home/" target="_blank">French Institute of Bioinformatics</a>`,
	description: "IFB Scenario 2",
	tags: ["unix", "shell", "terminal"],
	tools: ["ls", "date"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Manipulating files and directories", component: Step0 },
		{ name: "Tree, path & files", component: Step1 },
		{ name: "Absolute paths", component: Step2 },
		{ name: "Relative paths", component: Step3 },
		{ name: "Change directory", component: Step4 },
		{ name: "The HOME directory", component: Step5 },
		{ name: "Create or copy", component: Step6 },
		{ name: "Move or remove", component: Step7 },
		{ name: "Congratulations", component: Step8 }
	],
	init: "mkdir -p /root/tutorial/projects/facts",
	files: [
		"bos_taurus/UMD3.1/star-2.7.2b/chrStart.txt",
		"bos_taurus/UMD3.1/star-2.7.2b/chrLength.txt",
		"bos_taurus/UMD3.1/star-2.7.2b/SAindex",
		"bos_taurus/UMD3.1/star-2.7.2b/chrName.txt",
		"bos_taurus/UMD3.1/star-2.7.2b/genomeParameters.txt",
		"bos_taurus/UMD3.1/bowtie2/Bos_taurus.UMD3.1.dna.toplevel_F.2.bt2",
		"bos_taurus/UMD3.1/bowtie2/Bos_taurus.UMD3.1.dna.toplevel_F.1.bt2",
		"bos_taurus/UMD3.1/bowtie2/Bos_taurus.UMD3.1.dna.toplevel_F.rev.1.bt2",
		"bos_taurus/UMD3.1/bowtie2/Bos_taurus.UMD3.1.dna.toplevel_F.rev.2.bt2",
		"bos_taurus/UMD3.1/fasta/Bos_taurus.UMD3.1.dna.toplevel_F.fa",
		"bos_taurus/UMD3.1/fasta/Bos_taurus.UMD3.1.dna.toplevel_F.fa.fai",
		"nr/nr_2018-09-28/diamond/viral.protein_refseq_98.dmnd",
		"nr/nr_2018-09-28/diamond/nr.dmnd",
		"nr/nr_2018-09-28/fasta/nr.fsa",
		"nr/nr_2018-09-28/blast/nr.02.psd",
		"nr/nr_2018-09-28/blast/nr.01.psd",
		"nr/nr_2018-09-28/blast/nr.02.phd",
		"nr/nr_2018-09-28/blast/nr.01.phd",
		"nr/nr_2018-09-28/blast/nr.01.ppi",
		"nr/nr_2018-09-28/blast/nr.02.ppi",
		"homo_sapiens/hg19/bowtie2/hg19.rev.2.bt2",
		"homo_sapiens/hg19/bowtie2/hg19.2.bt2",
		"homo_sapiens/hg19/bowtie2/hg19.rev.1.bt2",
		"homo_sapiens/hg19/bowtie2/hg19.1.bt2",
		"homo_sapiens/hg19/hisat2/hg19.3.ht2",
		"homo_sapiens/hg19/hisat2/hg19.2.ht2",
		"homo_sapiens/hg19/hisat2/hg19.1.ht2",
		"homo_sapiens/hg19/hisat2/hg19.4.ht2",
		"homo_sapiens/hg19/fasta/hg19.fa.fai",
		"homo_sapiens/hg19/fasta/hg19.fa",
		"homo_sapiens/hg19/star-2.7.5a/chrStart.txt",
		"homo_sapiens/hg19/star-2.7.5a/chrLength.txt",
		"homo_sapiens/hg19/star-2.7.5a/SAindex",
		"homo_sapiens/hg19/star-2.7.5a/chrName.txt",
		"homo_sapiens/hg19/star-2.7.5a/genomeParameters.txt",
		"homo_sapiens/hg38/bowtie2/hg38.rev.2.bt2",
		"homo_sapiens/hg38/bowtie2/hg38.rev.1.bt2",
		"homo_sapiens/hg38/bowtie2/hg38.2.bt2",
		"homo_sapiens/hg38/bowtie2/hg38.1.bt2",
		"homo_sapiens/hg38/hisat2/genome.2.ht2",
		"homo_sapiens/hg38/hisat2/genome.3.ht2",
		"homo_sapiens/hg38/hisat2/genome.1.ht2",
		"homo_sapiens/hg38/hisat2/genome.4.ht2",
		"homo_sapiens/hg38/fasta/hg38.fa",
		"homo_sapiens/hg38/fasta/hg38.fa.fai",
		"homo_sapiens/hg38/star-2.7.5a/chrStart.txt",
		"homo_sapiens/hg38/star-2.7.5a/chrLength.txt",
		"homo_sapiens/hg38/star-2.7.5a/SAindex",
		"homo_sapiens/hg38/star-2.7.5a/chrName.txt",
		"homo_sapiens/hg38/star-2.7.5a/genomeParameters.txt",
		"test/first_file.txt",
		"test/second_file.txt"
	]
};

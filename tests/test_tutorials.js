// Test 1 representative command per tutorial

import { get } from "svelte/store";
import { CLI } from "../src/terminal/cli"; const $CLI = get(CLI);
let observed;

describe("Test tutorial contents (1 representative command)", () => {
	before(async () => {
		console.log("Initializing Aioli");
		await $CLI.init({
			tools: ["bedtools/2.29.2", "bowtie2/bowtie2-align-s/2.4.2", "samtools/1.10", "bcftools/1.10"],
			files: [
				// bedtools
				"public/data/bedtools-intro/exons.bed",
				"public/data/bedtools-intro/cpg.bed",
				"public/data/bedtools-intro/gwas.bed",
				"public/data/bedtools-intro/hesc.chromHmm.bed",
				// bowtie2
				"public/data/bowtie2-intro/reads_1.fq",
				"public/data/bowtie2-intro/reads_2.fq",
				// samtools
				"public/data/samtools-intro/sample.sam",
			]
		});
	});

	it("bedtools", async () => {
		observed = await $CLI.exec("bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted -wa -wb -names cpg gwas chromhmm | head -n 10000 | tail -n 10");
		expect(observed).to.equal("chr1\t935245\t935552\tNM_021170_exon_3_0_chr1_935246_r\t0\t-\tchromhmm\tchr1\t932537\t937537\t3_Poised_Promoter\nchr1\t948846\t948956\tNM_005101_exon_0_0_chr1_948847_f\t0\t+\tcpg\tchr1\t948670\t948894\tCpG:_19\nchr1\t948846\t948956\tNM_005101_exon_0_0_chr1_948847_f\t0\t+\tchromhmm\tchr1\t948337\t949337\t4_Strong_Enhancer\nchr1\t949363\t949919\tNM_005101_exon_1_0_chr1_949364_f\t0\t+\tcpg\tchr1\t949329\t949851\tCpG:_35\nchr1\t949363\t949919\tNM_005101_exon_1_0_chr1_949364_f\t0\t+\tchromhmm\tchr1\t949337\t949537\t2_Weak_Promoter\nchr1\t949363\t949919\tNM_005101_exon_1_0_chr1_949364_f\t0\t+\tchromhmm\tchr1\t949537\t949937\t6_Weak_Enhancer\nchr1\t955502\t955753\tNM_198576_exon_0_0_chr1_955503_f\t0\t+\tcpg\tchr1\t954768\t956343\tCpG:_148\nchr1\t955502\t955753\tNM_198576_exon_0_0_chr1_955503_f\t0\t+\tchromhmm\tchr1\t954537\t955537\t6_Weak_Enhancer\nchr1\t955502\t955753\tNM_198576_exon_0_0_chr1_955503_f\t0\t+\tchromhmm\tchr1\t955537\t956137\t2_Weak_Promoter\nchr1\t957580\t957842\tNM_198576_exon_1_0_chr1_957581_f\t0\t+\tchromhmm\tchr1\t957537\t958937\t9_Txn_Transition");
	});

	it("bowtie2", async () => {
		let stderr = "";
		observed = await $CLI.exec("REF=/bowtie2/example/index/lambda_virus; bowtie2 -x $REF -1 reads_1.fq -2 reads_2.fq -S eg2.sam", d => stderr = d);
		expect(observed).to.equal("");
		expect(stderr).to.equal("pthread_sigmask() is not supported: this is a no-op.\n25 reads; of these:\n  25 (100.00%) were paired; of these:\n    0 (0.00%) aligned concordantly 0 times\n    25 (100.00%) aligned concordantly exactly 1 time\n    0 (0.00%) aligned concordantly >1 times\n    ----\n    0 pairs aligned concordantly 0 times; of these:\n      0 (0.00%) aligned discordantly 1 time\n    ----\n    0 pairs aligned 0 times concordantly or discordantly; of these:\n      0 mates make up the pairs; of these:\n        0 (0.00%) aligned 0 times\n        0 (0.00%) aligned exactly 1 time\n        0 (0.00%) aligned >1 times\n100.00% overall alignment rate\n");
	});

	// Must run after bowtie2
	it("bcftools", async () => {
		let stderr = "";
		observed = await $CLI.exec("REF_FA=/bowtie2/example/reference/lambda_virus.fa; samtools view eg2.sam -o eg2.bam; samtools sort eg2.sam -o eg2.sorted.bam; bcftools mpileup -f $REF_FA eg2.sorted.bam | head -n2", d => stderr = d);
		expect(observed).to.equal(`##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">`);
		expect(stderr).to.equal(`[mpileup] 1 samples in 1 input files\n[mpileup] maximum number of reads per input file set to -d 250\n`);
	});

	it("samtools", async () => {
		observed = await $CLI.exec("samtools view -b sample.sam -o sample.bam; samtools sort sample.bam -o sample.sorted.bam; samtools index sample.sorted.bam; samtools view -c sample.sorted.bam 20:1.4e6-1.5e6");
		expect(observed).to.equal("338");
	});
});

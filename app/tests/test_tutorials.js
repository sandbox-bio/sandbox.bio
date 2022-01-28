// Test 1 representative command per tutorial

import { get } from "svelte/store";
import { CLI } from "../src/terminal/cli"; const $CLI = get(CLI);
import { TOOLS } from "./utils";

let observed;

describe("Test tutorial contents (1 representative command)", () => {
	before(async () => {
		console.log("Initializing Aioli");
		await $CLI.init({
			pwd: "tests",
			tools: ["bedtools/2.29.2", "bowtie2/bowtie2-align-s/2.4.2", "bcftools/1.10", "minimap2/2.22", "jq/1.6", "gawk/5.1.0", ...TOOLS],
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
				// gawk
				"public/data/awk-intro/orders.tsv",
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
		observed = await $CLI.exec("REF_FASTA=/bowtie2/example/reference/lambda_virus.fa; samtools view eg2.sam -o eg2.bam; samtools sort eg2.sam -o eg2.sorted.bam; bcftools mpileup -f $REF_FASTA eg2.sorted.bam", d => stderr = d);
		expect(observed).to.contain(`##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">`);
		expect(stderr).to.equal(`[mpileup] 1 samples in 1 input files\n[mpileup] maximum number of reads per input file set to -d 250\n`);
	});

	it("samtools", async () => {
		observed = await $CLI.exec("samtools view -b sample.sam -o sample.bam; samtools sort sample.bam -o sample.sorted.bam; samtools index sample.sorted.bam; samtools view -c sample.sorted.bam 20:1.4e6-1.5e6");
		expect(observed).to.equal("338");
	});

	it("minimap2", async () => {
		observed = await $CLI.exec("minimap2 -a /minimap2/MT-human.fa /minimap2/MT-orang.fa");
		expect(observed).to.contain("@SQ\tSN:MT_human\tLN:16569\n@PG\tID:minimap2\tPN:minimap2\tVN:2.22-r1101\tCL:minimap2 -a /minimap2/MT-human.fa /minimap2/MT-orang.fa\nMT_orang\t0\tMT_human\t577\t60\t14M2D4M3I37M1D85M1D232M1I559M1D6M1I550M1D2M1D146M2I3M1D3M1D132M1D3M1I40M3I13M1D1M1D335M3I4M1D3M2D342M1D52M1I13M3I1M2D52M1I592M1D3M1I485M1D5M1I974M3I4M3D230M1D59M1D156M1D31M1I98M1I26M14I329M3I7M3D1203M1D4M1I70M1D345M1D9M1I398M7I8M8I1M1I9M3I2M1I2M1D390M1I5M1D193M1I6M1D195M1D7M1I1826M1D10M1I1256M1D49M1D157M3D5M3I48M2I1M1I3M3D1203M1I2M2D1M1I44M2D2M1I2M1I38M2D16M2I2079M1I5M1D50M1D3M1I43M5I57M1I54M4D19M1I39M2D8M1I7M1I22M1I5M1I4M1D5M1I2M2D29M2I20M1D13M2I1M1D8M1D45M1D15M1I5M2D17M1D56M1D2M1I131M1I38M474S");
	});

	it("jq", async () => {
		observed = await $CLI.exec(`echo '{"test":{"something": "here"}}' | jq -r '.test.something'`);
		expect(observed).to.equal(`here`);

		// If don't specify `-r`, we get color output!
		observed = await $CLI.exec(`echo '{"test":{"something": "here"}}' | jq '.test.something'`);
		expect(observed).to.equal(`\u001b[0;32m"here"\u001b[0m`);
	});

	it("gawk", async () => {
		observed = await $CLI.exec(`awk -F "\t" ' { if($3 == "Chicken Bowl") sum += $2 } END { print(sum) }' orders.tsv`);
		expect(observed).to.equal(`761`);
	});
});

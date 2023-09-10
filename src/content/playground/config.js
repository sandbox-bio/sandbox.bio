export const config = {
	id: "playground",
	listed: false,
	tools: ["jq"],
	init: "clear",
	intro: `\u001b[0;37m# This playground is for open-ended exploration.\n# For guided tutorials, see https://sandbox.bio/tutorials\n#\n# Example:\n#   samtools view -o test.bam /shared/samtools/examples/toy.sam\n#   samtools index test.bam\n#   ls test.bam.bai\n#   samtools idxstats test.bam  # idxstats uses the .bai file \u001b[0m\n`,
	steps: []
};

# Exercise: Adapter Contamination Check

You've received a new batch of sequencing data and need to assess the level of Illumina adapter contamination.

<Exercise
id="ex1-adapter"
hints={[
"Use the -c flag to count matches",
"The Illumina TruSeq adapter sequence is AGATCGGAAGAGC",
"The file you need is reads.fastq"
]}
validate={result => {
return result.includes("40");
}}

>   <span slot="prompt">

    Count how many reads in `reads.fastq` contain the Illumina TruSeq adapter sequence `AGATCGGAAGAGC`.

  </span>
</Exercise>

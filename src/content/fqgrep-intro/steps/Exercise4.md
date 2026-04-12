# Exercise: Paired-End Adapter Screening

You need to find all read pairs where **either** R1 or R2 contains adapter contamination.

<Exercise
id="ex4-paired"
hints={[
"You need a flag that treats two input files as paired",
"The paired files are reads_R1.fastq and reads_R2.fastq"
]}
validate={result => {
return result.includes("55");
}}

>   <span slot="prompt">

    Count how many **read pairs** have Illumina adapter contamination (`AGATCGGAAGAGC`) in either R1 or R2. Use the files `reads_R1.fastq` and `reads_R2.fastq`.

  </span>
</Exercise>

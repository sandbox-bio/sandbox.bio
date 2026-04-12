# Exercise: Tn5 Mosaic End Detection

ATAC-seq and many tagmentation-based library preps leave Tn5 transposase mosaic end sequences in the reads.

<Exercise
id="ex2-tn5"
hints={[
"The Tn5 mosaic end sequence is AGATGTGTATAAGAGACAG",
"Use fqgrep with the -c flag",
"Search in reads.fastq"
]}
validate={result => {
return result.includes("35");
}}

>   <span slot="prompt">

    Count how many reads contain the Tn5 mosaic end sequence `AGATGTGTATAAGAGACAG`.

  </span>
</Exercise>

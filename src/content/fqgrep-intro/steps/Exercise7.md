# Exercise: IUPAC Degenerate Search

Your colleague designed a primer with the degenerate sequence `ACGTNNNNTGCA` and wants to know how many reads contain it. Rather than translating the IUPAC codes to regex by hand, use fqgrep's built-in IUPAC support.

<Exercise
id="ex7-iupac"
hints={[
"The --iupac flag requires a specific matching mode to be enabled first",
"Check the IUPAC step for the available modes"
]}
validate={result => {
return result.includes("20");
}}

>   <span slot="prompt">

    Use fqgrep's `--iupac` flag to count reads in `reads.fastq` matching the degenerate barcode `ACGTNNNNTGCA`.

  </span>
</Exercise>

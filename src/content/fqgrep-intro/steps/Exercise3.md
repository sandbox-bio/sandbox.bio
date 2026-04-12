# Exercise: Degenerate Barcode Pattern

Your library uses barcodes with the structure `ACGT-NNNN-TGCA` where `NNNN` represents 4 variable bases (a "degenerate" region).

<Exercise
id="ex3-regex"
hints={[
"You need a regular expression to match the variable region",
"Remember which regex character matches any single character"
]}
validate={result => {
return result.includes("20");
}}

>   <span slot="prompt">

    Count how many reads in `reads.fastq` contain a barcode matching `ACGT` followed by any 4 bases followed by `TGCA`.

  </span>
</Exercise>

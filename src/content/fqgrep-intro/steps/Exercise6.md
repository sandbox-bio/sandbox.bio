# Exercise: Protein Motif Search

You have a FASTQ file from a protein sequencing experiment and need to check for polyhistidine (His) tag contamination.

<Exercise
id="ex6-protein"
hints={[
"fqgrep needs to know the input contains amino acids, not DNA",
"The protein data file is proteins.fastq"
]}
validate={result => {
return result.includes("8");
}}

>   <span slot="prompt">

    Count how many reads in `proteins.fastq` contain a 6xHis tag (`HHHHHH`).

  </span>
</Exercise>

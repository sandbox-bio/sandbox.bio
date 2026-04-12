# Exercise: Thorough Contaminant Screen

Your lab's `contaminants.txt` file lists contaminant sequences in forward orientation only. But depending on library prep, adapters and Tn5 sequences can appear in either strand.

Perform a thorough screen that accounts for **both orientations** of all contaminants, and count the reads that pass.

<Exercise
id="ex8-thorough"
hints={[
"You need to combine three concepts: inverted matching, reverse complement, and a pattern file",
"Think about which flags from earlier steps handle each part"
]}
validate={result => {
return result.includes("417");
}}

>   <span slot="prompt">

    How many reads in `reads.fastq` are free of all contaminant sequences from `contaminants.txt` when checking **both** forward and reverse complement orientations?

  </span>
</Exercise>

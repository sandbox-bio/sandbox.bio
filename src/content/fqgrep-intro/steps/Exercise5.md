# Exercise: Extract Reads by Name

A colleague has given you a list of read names in `read_names.txt` that they want to inspect more closely. Extract just those reads and count them.

<Exercise
id="ex5-names"
hints={[
"There's a flag specifically for filtering by read name from a file",
"Check the Read Name Filtering step for the right flag"
]}
validate={result => {
return result.includes("5");
}}

>   <span slot="prompt">

    Count how many reads in `reads.fastq` match the names listed in `read_names.txt`.

  </span>
</Exercise>

<script>
import Execute from "$components/Execute.svelte";
</script>

We can now update our script file from previous tutorials:

<Execute command="nano getproteins.sh" />

To contain the following lines:

```bash
# CHEBI identifier given as input is renamed to ID
ID=$1

# Removes any previous files
rm -f chebi_${ID}_*.xml

csv_file=chebi_${ID}_xrefs_UniProt_relevant_identifiers.csv

cat "$csv_file" | \\
  xargs \\
    -I {} curl \\
    -o chebi_${ID}_{}.xml \\
    'https://rest.uniprot.org/uniprotkb/{}.xml'
```

#### Variable

We should note that the last line now includes the `xargs` and `curl` commands, and the `$ID` variable. This new variable is created in the first line to
contain the first value given as argument (`$1`). So, every time we mention
`$ID` in the script we are mentioning the first value given as argument. This
avoids ambiguity in cases where `$1` is used for other purposes. Since the preceding
character of `$ID` is an underscore (`_`), we have to add a backslash (`\`)
before it. The second line uses the rm command to remove any files that were
downloaded in a previous execution. We also now added two comments after
`#` character, so we humans do not forget why these commands are needed
for.

To execute the script once more:

<Execute command="chmod u+x getproteins.sh" />

<Execute command="./getproteins.sh 27732" />

And again, to check the results:

<Execute command="head -n 1 chebi_27732_*.xml | less" />

Type `q` to exit from `less`.

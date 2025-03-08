<script>
import Execute from "$components/Execute.svelte";
</script>

#### Data elements selection

Now we need to select just the first column, the one that contains the protein
identifiers. Selecting columns from a tabular file is one easy task for `cut`.
The cut command can receive as arguments the character that divides each
data element (column) in a line using the `-d` option, and the `-f` option to
indicate which columns to select. The equivalent long form to the `-d` option
is `--delimiter=DELIM`. The equivalent long form to the `-f` option is `--fields=LIST`.
For example, we can get the first column of our CSV file:

<Execute command="cut -d, -f1 < chebi_27732_xrefs_UniProt_relevant.csv" />

We should note that comma (`,`) is the character that separates data elements
in a CSV file, and represents the first data element.

The command will display only the first column of the file, i.e. the protein
identifiers.

For example, we can get the first and third columns separated by a comma:

<Execute command="cut -d, -f1,3 < chebi_27732_xrefs_UniProt_relevant.csv" />

Now, the output contains both the first and third column of the file.

We can update our script file:

<Execute command="nano getproteins.sh" />

To contain the following lines:

```bash
URL="https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-e=1&6578706f7274=1&chebiId=$1&dbName=UniProt"

curl -s "$URL" | \\
    grep -e 'CC - MISCELLANEOUS' \\
        -e 'CC - DISRUPTION PHENOTYPE' \\
        -e 'CC - DISEASE' | \\
    cut -d, -f1
```

The last line is the only that changes, except the `| \` in the previous line to
redirect the output.
To execute the script, we can type again:

<Execute command="./getproteins.sh 27732" />

The output should be similar of what we got previously, but now only the
protein identifiers are displayed.

To save the output as a file with the relevant proteins’ identifiers, we only
need to add the redirection operator:

<Execute command="./getproteins.sh 27732 > chebi_27732_xrefs_UniProt_relevant_identifiers.csv" />

And we may now check the contents of the created file:

<Execute command="cat chebi_27732_xrefs_UniProt_relevant_identifiers.csv" />

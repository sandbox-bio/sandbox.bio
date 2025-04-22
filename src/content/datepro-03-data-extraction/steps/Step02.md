<script>
import Execute from "$components/Execute.svelte";
</script>

To use multiple patterns, we must precede each pattern with the `-e` option:

<Execute command="grep -e 'CC - MISCELLANEOUS' -e 'CC - DISRUPTION PHENOTYPE' -e 'CC - DISEASE' chebi_27732_xrefs_UniProt.csv " />

The equivalent long form to the `-e` option is `--regexp=PATTERN`.

The output on our terminal should be a longer list of proteins.

We should note that as previously, we can add `| less` to check all of
them more carefully. The less command also gives the opportunity to find
lines based on a pattern. We only need to type `/` and then a pattern. To exit from `less`, type `q`.

<Execute command="grep -e 'CC - MISCELLANEOUS' -e 'CC - DISRUPTION PHENOTYPE' -e 'CC - DISEASE' chebi_27732_xrefs_UniProt.csv | less" />

We can now update our script file

<Execute command="nano getproteins.sh" />

To contain the following lines:

```bash
url="https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-e=1&6578706f7274=1&chebiId=$1&dbName=UniProt"

curl -s "$url" | \
    grep \
      -e 'CC - MISCELLANEOUS' \
      -e 'CC - DISRUPTION PHENOTYPE' \
      -e 'CC - DISEASE'
```

> The `curl` command here at this platform will be replaced by a local version to access the data locally instead of online due to access restrictions. However, the same command will work just fine in your local terminal.

We should note that we added the `-s` option to suppress the progress information of curl, and the characters `| \` to the end of line to redirect the
output of that line as input of the next line, in this case the grep command.
We need to be careful in ensuring that `\` is the last character in the line, i.e.
spaces in the end of the line may cause problems.

Add the right permissions with `chmod` as we did previously:

<Execute command="chmod u+x getproteins.sh" />

We can now execute the script again:

<Execute command="./getproteins.sh 27732" />

The output should be similar of what we got previously, but the script downloads the data and filters immediately. To save the file with the relevant proteins, we only need to add the redirection operator:

<Execute command="./getproteins.sh 27732 > chebi_27732_xrefs_UniProt_relevant.csv" />

The next step will show how to extract only the protein identifiers from the previous list.

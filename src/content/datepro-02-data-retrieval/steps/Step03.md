<script>
import Execute from "$components/Execute.svelte";
</script>

If we want to analyze all the lines we can redirect the output to the command line tool less, which allows us to navigate through the output by using the arrow keys. To do that we can add the bar character (`|`) between two commands, which will transfer the output of the first command as input of the second:

<Execute command="./getproteins.sh 27732 | less" />

To exit from less just press **q**.

However, what we really want is to save the output as a file, not just
printing some characters on the screen. Thus, what we should do is redirect
the output to a CSV file. This can be done by adding the redirect operator `>`
and the filename, as described previously:

<Execute command="./getproteins.sh 27732 > chebi_27732_xrefs_UniProt.csv" />

We should note that `curl` still prints some progress information into the
terminal.

<Alert>In your local term the progress information will more detailed.</Alert>

#### Standard error output

This happens because it is displaying that information into the standard er-
ror output, which was not redirected to the file. The `>` character without
any preceding number by default redirects the standard output. The same
happens if we precede it by the number 1. If we do not want to see that
information, we can also redirect the standard error output (2), but in this
case to the null device (/dev/null):

<Execute command="./getproteins.sh 27732 > chebi_27732_xrefs_UniProt.csv 2>/dev/null" />

We can also use the -s option of curl in order to suppress the progress
information,:

```bash
curl -s "https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-=1&6578706f7274=1&chebiId=$1&dbName=UniProt"
```

by adding it to our script:

<Execute command="nano getproteins.sh" />


Now when executing the script, no progress information is shown:

<Execute command="./getproteins.sh 27732 > chebi_27732_xrefs_UniProt.csv" />

To check if the file was really created and to analyze its contents, we can
use the less command:

<Execute command="less chebi_27732_xrefs_UniProt.csv" />

We can also open the file in our spreadsheet application, such as LibreOffice
Calc or Microsoft Excel.

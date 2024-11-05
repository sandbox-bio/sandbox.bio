<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

Now that we have the URLs, we can automatically download the files using
the curl command instead of echo:

<Execute command="cat chebi_27732_xrefs_UniProt_relevant_identifiers.csv | xargs -I &lcub;&rcub; curl -o 'chebi_27732_&lcub;&rcub;.xml' 'https://rest.uniprot.org/uniprotkb/&lcub;&rcub;.xml' " />

<Alert>The `curl` command on this platform has access restrictions but works with `rest.uniprot.org` links.</Alert>

We should note that we now use the `-o` option to save the output to a given
file, named after each protein identifier. The equivalent long form to the `-o`
option is `--output <file>`.
To check if everything worked as expected we can use the ls command to
view which files were created:

<Execute command="ls chebi_27732_*.xml" />

The asterisk character (`*`) is here used to represent any file whose name
starts with `chebi_27732_` and ends with `.xml`.
To check the contents of any of them, we can use the less command:

<Execute command="less chebi_27732_P21817.xml" />

#### File header

We should note that the content of every file has to start with `<?xml` otherwise there was a download error, and we have to run curl again for those
entries. To check the header of each file, we can use the head command
together with `less`.

<Execute command="head -n 1 chebi_27732_*.xml | less" />

To exit from `less`, type `q`.

The `-n` option specifies how many lines to print, in the previous command
just one.
If for any reason, we are not able to download the files from UniProt, we
can get them from the [book file archive](http://labs.rd.ciencias.ulisboa.pt/book/)

In the next step, we will create a script using `xargs`.

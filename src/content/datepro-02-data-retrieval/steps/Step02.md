<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

#### Data Retrieval

After having the link, we need a web retrieval tool that works like our internet
browser, i.e. receives as input a URL for programmatic access and retrieves
its contents from the internet. We will use Client Uniform Resource Locator
(cURL), which is available as a command line tool, and allows us to download
the result of opening a URL directly into a file.

<Alert>The `curl` command here at this platform will be replaced by a local version to access the data locally instead of online due to access restrictions. However, the same command will work just fine in your local terminal.</Alert>

For example, to display in our screen the list of proteins related to caffeine,
we just need to add the respective URL as input argument:

<Execute command="curl 'https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-e=1&6578706f7274=1&chebiId=27732&dbName=UniProt'" />

The output on our terminal is the long list of proteins.

An alternative to curl is the command `wget` (not available at this plaform), which also receives a URL as argument but by default `wget` writes the contents to a file instead of displaying it on the screen. So, the equivalent command, is to add the -O- option to select where the contents is placed.

```bash
wget -O- 'https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-=1&6578706f7274=1&chebiId=27732&dbName=UniProt'
```

Instead of using a fixed URL, we can update the script:

<Execute command="nano getproteins.sh" />

To contain only the following line:

```bash
curl "https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-=1&6578706f7274=1&chebiId=$1&dbName=UniProt"
```

We should note that now we are using double quotes, since we replaced the
caffeine identifier by `$1`.

Now to execute the script we only need to provide a ChEBI identifier as input argument:

<Execute command="./getproteins.sh 27732" />

The output on our terminal is the long list of proteins.

Or, if we want the proteins related to carbon monoxide, we only need to
replace the argument:

<Execute command="./getproteins.sh 17245" />

And the output on our terminal should be an even longer list of proteins.

The next step will show how to manage these long lists.

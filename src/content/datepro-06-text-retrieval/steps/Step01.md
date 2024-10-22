<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

Now that we have all the PubMed identifiers, we need to download the text
included in the titles and abstracts of each publication.

#### Publication URL

To retrieve from the UniProt citations service the publication entry of a given
identifier, we can again use the `curl` command and a link to the publication
entry. For example, if we click on the Format button of the UniProt citations
service entry, we can get the link to the RDF/XML version. RDF is a
standard data model that can be serialized in a XML format. Thus, in our
case, we can deal with this format like we did with XML.
We can retrieve the publication entry by executing the following command:

<Execute command="curl https://rest.uniprot.org/citations/1354642.rdf" />

<Alert>The `curl` command on this platform has access restrictions but works with  `rest.uniprot.org` and `eutils.ncbi.nlm.nih.gov` links.</Alert>

Alternatively, we can use the web service provided by PubMed at NCBI ,
by still using curl but with another link:

<Execute command="curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=1354642&retmode=text&rettype=xml'" />

The result is in XML and we can replace the PubMed identifier `135464` by a
comma separated list of identifiers, such has `2298749,1354642,8220422`.
Thus, we can now update the script 

<Execute command="nano getpublications.sh" />

To have the following commands:

<pre class="code border p-2" style="white-space: pre-wrap">
ID=$1 # The CHEBI identifier given as input is renamed to ID
rm -f chebi\_$&lcub;ID&rcub;\_*.rdf # Removes any previous files

grep -l '&lt;name type="scientific">Homo sapiens&lt;/name>' chebi\_$&lcub;ID&rcub;\_*.xml | \\
    xargs -I &lcub;&rcub; grep '&lt;dbReference type="PubMed"' &lcub;&rcub; | \\
    cut -d'"' -f4 | sort -u | \\
    xargs -I &lcub;&rcub; curl -o chebi\_$&lcub;ID&rcub;\_&lcub;&rcub;.rdf 'https://rest.uniprot.org/citations/&lcub;&rcub;.rdf' 
</pre>

Again, do not forget to save it in our working directory, and add the right
permissions with chmod as we did previously with the other scripts.

<Execute command="chmod u+x getpublications.sh" />

To execute the script again:

<Execute command="./getpublications.sh 27732" />

It may take a while to download all the entries, but probably no more than
one minute with a standard internet connection.

To check if everything worked as expected we can use the ls command to
view which files were created:

<Execute command="ls chebi_27732_*.rdf" />

If for any reason, we are not able to download the abstracts from UniProt,
we can get them from the [book file archive](http://labs.rd.ciencias.ulisboa.pt/book/).

In the next step, we will extract the titles and abstracts from the RDF files.

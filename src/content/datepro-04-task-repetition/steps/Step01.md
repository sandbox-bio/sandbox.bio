<script>
import Execute from "$components/Execute.svelte";
</script>

#### Task Repetition

Given a protein identifier we can construct the URL that will enable us to
download its information from UniProt. We can use the RESTful web services
provided by UniProt 39 , more specifically the one that allow us to retrieve a
specific entry 40 . The construction of the URL is simple, it starts always
by https://rest.uniprot.org/uniprotkb/, followed by the protein
identifier, ending with a dot and the data format. For example, the link for
protein P21817 using the XML format is:
https://rest.uniprot.org/uniprotkb/P21817.xml

###### Assembly line

However, we need to construct one URL for each protein from the list we
previously retrieved. The size of the list can be large (hundreds of proteins),
varies for different compounds and evolves with time. Thus, we need an
assembly line in which a list of proteins identifiers, independently of its size,
are added as input to commands that construct one URL for each protein and
retrieve the respective file.

The `xargs` command line tool works as an assembly line, it executes a
command per each line given as input.
We can start by experimenting the `xargs` command by giving as input the list of protein identifiers in file `chebi_27732_xrefs_UniProt_relevant_identifiers.csv`,
and display each identifier on the screen in the middle of a text message by
providing the echo command as argument:

<Execute command="cat chebi_27732_xrefs_UniProt_relevant_identifiers.csv | xargs -I {} echo 'Another protein id {} to retrieve'" />

The xargs command received as input the contents our CSV file, and for
each line displayed a message including the identifier in that line. The `-I`
option tells `xargs` to replace `{}` in the command line given as argument by
the value of the line being processed. The equivalent long form to the `-I`
option is `--replace=R`.

Instead of creating inconsequential text messages, we can use xargs to
create the URLs:

<Execute command="cat chebi_27732_xrefs_UniProt_relevant_identifiers.csv | xargs -I {} echo 'https://rest.uniprot.org/uniprotkb/{}.xml'" />

We can try to use these links in our internet browser to check if those
displayed URLs are working correctly.

In the next step, we will automatically download the files using the `curl` command.

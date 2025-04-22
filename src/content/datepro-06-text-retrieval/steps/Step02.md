<script>
import Execute from "$components/Execute.svelte";
</script>

#### Title and Abstract

Each file has the title and abstract of the publication as values of the title
and `rdfs:comment` elements, respectively. To extract them we can again
use the `xmllint` command:

<Execute command={`xmllint --xpath "//*[local-name()='title' or local-name()='comment']" chebi_27732_1354642.rdf`} />

The output should be the text inside XML elements.
To remove the XML elements, we can again add `text()` to the XPath
query:

<Execute command={`xmllint --xpath "//*[local-name()='title' or local-name()='comment']/text()" chebi_27732_1354642.rdf`} />

The output should now be free of XML elements.

Thus, let us create the script `gettext.sh`

<Execute command="nano gettext.sh" />

To have the following commands:

```bash
# CHEBI identifier given as input is renamed to ID
ID=$1

xmllint \
    --xpath "//*[local-name()='title' or local-name()='comment']/text()" chebi_$ID_*.rdf
```

Again do not forget to save it in our working directory, and add the right
permissions.

<Execute command="chmod u+x gettext.sh" />

Now to execute the script and see the retrieved text:

<Execute command="./gettext.sh 27732 | less" />

Type `q` to exit from `less`.

We can save the resulting text in a file named `chebi_27732.txt` that we
may share or read using our favorite text editor, by adding the redirection
operator:

<Execute command="./gettext.sh 27732 > chebi_27732.txt" />

In the next step, we will recognize diseases in the extracted text.

<script>
import Execute from "$components/Execute.svelte";
</script>

Assuming that our information need only concerns human diseases, we have
to process the XML file of each protein to check if it represents a Homo sapiens
(Human) protein.
Human proteins
For performing this filter, we can again use the grep command, to select only
the lines of any XML file that specify the organism as Homo sapiens:

<Execute command="grep '<name type="scientific">Homo sapiens</name>' chebi_27732_*.xml" />

We should get in our display the filenames that represent a human protein.

We should note that since the asterisk character (`*`) provides multiple files
as argument to grep, the ones whose name starts with `chebi_27732_` and
ends with `.xml`, the output now includes the filename (followed by a colon)
where each line was matched.
We can use the `cut` command to extract only the filename, but grep has
the `-l` option to just print the filename:

<Execute command="grep -l '<name type="scientific">Homo sapiens</name>' chebi_27732_*.xml" />

The equivalent long form to the `-l` option is `--files-with-matches`.

The output will now show only the filenames. These four files represent the four Human proteins related to caffeine.

##### PubMed identifiers

Now we need to extract the PubMed identifiers from these files to retrieve
the related publications. For example, if we execute the following command:

<Execute command="grep 'dbReference type="PubMed"' chebi_27732_P21817.xml" />

The output is a long list of publications related to protein P21817.

In the next step, we will learn how to extract these identifiers from all human proteins. 
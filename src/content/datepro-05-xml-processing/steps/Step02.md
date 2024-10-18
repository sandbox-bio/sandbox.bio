<script>
import Execute from "$components/Execute.svelte";
</script>

To extract just the identifier, we can again use the cut command:

<Execute command={`grep 'dbReference type="PubMed"' chebi_27732_P21817.xml | cut -d'"' -f4`} />

We should note that `"` is used as the separation character and, since the
PubMed identifier appears after the third `"`, the `4` represents the identifier.
Now the output should be just the list of identifiers.

<Execute command={`grep '&lt;name type="scientific">Homo sapiens&lt;/name>' chebi_27732_*.xml`} />

#### PubMed identifiers extraction

Now to apply to every protein we may again use the `xargs` command:

<Execute command={`grep -l '&lt;name type="scientific">Homo sapiens&lt;/name>' chebi_27732_*.xml | \\ xargs -i &lcub;&rcub; grep 'dbReference type="PubMed"' &lcub;&rcub; | \\ cut -d'"' -f4`} />

This may provide a long list of PubMed identifiers, including repetitions since
the same publication can be cited in different entries.

#### Duplicate removal

To help us identify the repetitions, we can add the sort command, which will display the repeated identifiers in consecutive lines (due by sorting all identifiers):

<Execute command={`grep -l '&lt;name type="scientific">Homo sapiens&lt;/name>' chebi_27732_*.xml | \\ xargs -i &lcub;&rcub; grep 'dbReference type="PubMed"' &lcub;&rcub; | \\ cut -d'"' -f4 | \\ sort`} />

For example some repeated PubMed identifiers that we should easily be
able to see, such as `9607712`.

Fortunately, we also have the `-u` option that removes all these duplicates:

<Execute command={`grep -l 'name type="scientific">Homo sapiens&lt;/name>' chebi_27732_*.xml | xargs -i &lcub;&rcub; grep 'dbReference type="PubMed"' &lcub;&rcub; | cut -d'"' -f4 | sort -u`} />

To easily check how many duplicates were removed, we can use the word
count `wc` command with and without the usage of the `-u` option:

<Execute command={`grep -l 'name type="scientific">Homo sapiens&lt;/name>' chebi_27732_*.xml | xargs -i &lcub;&rcub; grep 'dbReference type="PubMed"' &lcub;&rcub; | cut -d'"' -f4 | sort | wc`} />

<Execute command={`grep -l 'name type="scientific">Homo sapiens&lt;/name>' chebi_27732_*.xml | xargs -i &lcub;&rcub; grep 'dbReference type="PubMed"' &lcub;&rcub; | cut -d'"' -f4 | sort -u | wc`} />

In case we have in our folder any auxiliary file, such as `chebi_27732_P21817_entry.xml`, we should add the option --exclude *entry.xml to the first `grep` command.

`wc` prints the numbers of lines, words, and bytes, thus in our case we are
interested in first number.
We can see that we have removed 263 − 133 = 130 duplicates.
Just for curiosity, we can also use the shell to perform simple mathematical
calculations using the `expr` command:

<Execute command="expr 263 - 133" />

Now let us create a script file named getpublications.sh by using a text
editor:

<Execute command="nano getpublications.sh" />

to add the following lines:

<pre class="code border p-2" style="white-space: pre-wrap">ID=$1 # The CHEBI identifier given as input is renamed to ID
grep -l '&lt;name type="scientific">Homo sapiens&lt;/name>' chebi\_$ID\_*.xml | \
xargs -i &lcub;&rcub; grep '&lt;dbReference type="PubMed"' &lcub;&rcub; | \
cut -d'"' -f4 | sort -u
</pre>

Again, do not forget to save it in our working directory, and add the right
permissions with chmod as we did previously with the other scripts.

<Execute command="chmod u+x getpublications.sh" />

To execute the script again:

<Execute command="./getpublications.sh 27732" />

We can verify how many unique publications were obtained by using the
`-l` option of `wc`, that provides only the number of lines:

<Execute command="./getpublications.sh 27732 | wc -l"  />

The output will be 133 as expected.

In the next step, we will learn how to perform more complex queries in XML files. 
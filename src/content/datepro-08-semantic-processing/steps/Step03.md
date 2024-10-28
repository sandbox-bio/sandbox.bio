<script>
import Execute from "$components/Execute.svelte";
</script>

## Synonyms

For example, to find all the synonyms of a disease, we can use the same XPath as used before but replacing the keyword label by `hasExactSynonym`:

<Execute command={`xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/DOID_8545']/*[local-name()='hasExactSynonym']" doid.owl`} />

The output will be the two synonyms of malignant hyperthermia

We can also get both the primary label and the synonyms. We only need to add an alternative match to the keyword label:

<Execute command={`xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/DOID_8545']/*[local-name()='hasExactSynonym' or local-name()='label']" doid.owl`} />

The output will include now the two synonyms plus the official label

Thus, we can now update the script `getlabels.sh`:

<Execute command="nano getlabels.sh" />

To include synonyms:

<pre class="code border p-2" style="white-space: pre-wrap">
OWLFILE=$1
xargs -I &lcub;&rcub; xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='&lcub;&rcub;']/*[local-name()='hasExactSynonym' or local-name()='hasRelatedSynonym' or local-name()='label']/text()" $OWLFILE</pre>

We can test the script exactly in the same way as before:

<Execute command="echo -e 'http://purl.obolibrary.org/obo/DOID_8545' | ./getlabels.sh doid.owl" />

But now the output will display multiple labels for this class.

#### URI of synonyms

Since the script now returns alternative labels, we may encounter some problems if we send the output to the `geturi.sh` script:

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_8545' | ./getlabels.sh doid.owl | ./geturi.sh doid.owl" />

The previous command will display XPath warnings for the two synonyms. If we do not want to know about these mismatches, we can always redirect them to the null device:

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_8545' | ./getlabels.sh doid.owl | ./geturi.sh doid.owl 2>/dev/null" />

However, we can update the script `geturi.sh`:
<Execute command="nano geturi.sh" />

To also include synonyms:

<pre class="code border p-2" style="white-space: pre-wrap">
OWLFILE=$1
xargs -I &lcub;&rcub; xmllint --xpath "//*[(local-name()='hasExactSynonym' or local-name()='hasRelatedSynonym' or local-name()='label') and text()='&lcub;&rcub;']/../@*[local-name()='about']" $OWLFILE | \
cut -d\" -f2</pre>

Now we can execute the same command:

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_8545' | ./getlabels.sh doid.owl | ./geturi.sh doid.owl" />

Every label should now be matched exactly with the same class.

If we want to avoid repetitions, we can add the sort command with the `-u` option to the end of each command, as we did previously.

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_8545' | ./getlabels.sh doid.owl | ./geturi.sh doid.owl | sort -u" />

The output should now be only one URI.

<script>
import Execute from "$components/Execute.svelte";
</script>

## My Lexicon

Now that we know how to extract all the labels and related classes from an ontology, we can construct our own lexicon with the list of terms that we want to recognize in text.
Let us start by creating the file `do_8545_lexicon.txt` representing our lexicon for malignant hyperthermia with all its labels:

<Execute command="echo 'malignant hyperthermia' | ./geturi.sh doid.owl | ./getlabels.sh doid.owl > do_8545_lexicon.txt" />

#### Ancestors labels

Now we can add to the lexicon all the labels of the ancestors of malignant hyperthermia by adding the redirection operator:

<Execute command="echo 'malignant hyperthermia' | ./geturi.sh doid.owl | ./getancestors.sh doid.owl | ./getlabels.sh doid.owl >> do_8545_lexicon.txt" />

We should note that now we use `>>` and not `>`, this will append more lines to the file instead of creating a new file from scratch.
Now we can check the contents of the file `do_8545_lexicon.txt` to see the terms we got:

<Execute command="cat do_8545_lexicon.txt | sort -u" />

We should note that we use the sort command with the `-u` option to eliminate any duplicates that may exist.

We can also apply the same commands for caffeine to produce its lexicon in the file `chebi_27732_lexicon.txt` by adding the redirection operator:

<Execute command="echo 'caffeine' | ./geturi.sh chebi_lite.owl | ./getlabels.sh chebi_lite.owl > chebi_27732_lexicon.txt" />

<Execute command="echo 'caffeine' | ./geturi.sh chebi_lite.owl | ./getancestors.sh chebi_lite.owl | ./getlabels.sh chebi_lite.owl >> chebi_27732_lexicon.txt" />

> Please note that it may take some time to retrieve all labels, even with these reduced OWL files. In the meantime, you can look at the following steps while it is running.

Now let us check the contents of this new lexicon:

<Execute command="cat chebi_27732_lexicon.txt | sort -u" />

Now we should be able to see that this lexicon is much larger.

#### Merging labels

If we are interested in finding everything related to caffeine or malignant hyperthermia, we may be interested in merging the two lexicons in a file named `lexicon.txt`:

<Execute command="cat do_8545_lexicon.txt chebi_27732_lexicon.txt | sort -u > lexicon.txt" />

Using this new lexicon, we can recognize any mention in our previous file
named `chebi_27732_sentences.txt`:

<Execute command="grep -w -i -F -f lexicon.txt chebi_27732_sentences.txt" />

We added the `-F` option because our lexicon is a list of fixed strings, i.e. does not include regular expressions. The equivalent long form to the `-F` option is `--fixed-strings`.
We now get more sentences, including some that do not include a direct
mention to caffeine or malignant hyperthermia. For example, the following sentence was selected because it mentions molecule, which is an ancestor of caffeine.

```text
The remainder of the molecule is hydrophilic
and presumably constitutes the cytoplasmic
domain of the protein.
```

Another example is the following sentence, which was selected because it mentions disease, which is an ancestor of malignant hyperthermia:

```text
Our data suggest that divergent activity profiles
may cause varied disease phenotypes by specific mutations.
```

We can also use our script `getentities.sh` (previous tutorial) giving this lexicon as argument. However, since we are not using any regular expressions it would be better to replace the `-E` option by `-F` to the `grep` command in the script, so the lexicon is interpreted as list of fixed strings to be
matched. Only then we can execute the script safely.

#### Ancestors matched

Besides these two previous examples, we can check if there other ancestors being matched by using the grep command with the `-o` option:

<Execute command="grep -o -w -F -f lexicon.txt chebi_27732_sentences.txt | sort -u" />

We can see that besides the terms caffeine and malignant hyperthermia,
only one ancestor of each one of them was matched, molecule and disease, respectively.

This can be explained because our text is somehow limited and because
we are using the official labels and we may be missing acronyms, and simple variations such as the plural of a term. To cope with this issue, we may use a [stemmer](https://en.wikipedia.org/wiki/Stemming), or use all the ancestors besides subsumption. However, if our lexicon is small is better to do it manually and maybe add some regular expressions to deal with some of the variations.

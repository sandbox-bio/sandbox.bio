<script>
import Execute from "$components/Execute.svelte";
</script>

## Entity Linking

When we are using a generic lexicon, we may be interested in identifying
what the recognized labels represent. For example, we may not be aware of what the matched label AD2 represents.

To solve this issue, we can use our script `geturi.sh` to perform entity linking (aka entity disambiguation, entity mapping, normalization), i.e. find the classes in the disease ontology that may be represented by the recognized label. For example, to find what AD2 represents, we can execute the following command:

<Execute command="echo 'AD2' | ./geturi.sh doid.owl" />

Only one URI is displayed.

Now we can retrieve other labels:

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_0110035' | ./getlabels.sh doid.owl" />

In this case, the result clearly shows that AD2 represents the _Alzheimer
disease_.

#### Modified labels

However, we may not be so lucky with the labels that were modified by our previous fixes in the lexicon. For example, we can test the case of ATR:

<Execute command="echo 'ATR' | ./geturi.sh doid.owl" />

As expected, we received the warning that no URI was found.

```text
XPath set is empty
```

An approach to address this issue may involve keeping a track of the original label in a lexicon using another file.

#### Ambiguity

We may also have to deal with ambiguity problems where a label may represent multiple terms. For example, if we check how many classes the acronym KOS may represent:

<Execute command="echo 'KOS' | ./geturi.sh doid.owl" />

We can see that it may represent two classes.

These two classes represent two distinct diseases, namely Kaufman oculocerebrofacial syndrome (DOID:0111456) and Kagami-Ogata syndrome (DOID:0111712), respectively.

We can also obtain their alternative labels by providing the two URI as
standard input to the `getlabels.sh` script:

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_0111456' | ./getlabels.sh doid.owl" />

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_0111712' | ./getlabels.sh doid.owl" />

We will get the following two lists, both containing KOS as expected.

If we find a `KOS` mention in the text, the challenge is to identify which of the syndromes the mention refers to. For addressing this challenge, we may have to use advanced entity linking techniques that analyze the context of the text.

#### Surrounding entities

An intuitive solution is to select the class closer in terms of meaning to the other classes mentioned in the surrounding text. This assumes that entities present in a piece of text are somehow semantically related to each other, which is normally the case. At least the author assumed some type of relation between them, otherwise the entities would not be in the same sentence.

Let us consider the following sentence about KOS:

```text
KOS is a syndromic intellectual disability
```

To identify the diseases in the previous sentence, we can execute the fol
lowing command:

<Execute command="echo 'KOS is a syndromic intellectual disability' | grep -o -w -F -f diseases.txt" />

We have a list of labels that can help us decide which is the right class
representing KOS.

To find their URIs we can use the `geturi.sh` script:

<Execute command="echo 'KOS is a syndromic intellectual disability' | grep -o -w -F -f diseases.txt | ./geturi.sh doid.owl" />

The only ambiguity is for KOS that returns two URIs, one representing the _Kaufman oculocerebrofacial syndrome _(DOID:0111456) and the other representing the _Kagami-Ogata syndrome_ (DOID:0111712).
The other URI represents the _Syndromic intellectual disability_ (DOID:0050888).

To decide which of the two URIs we should select, we can measure how
close in meaning they are to the other diseases also found in the text.

#### Semantic similarity

Semantic similarity measures have been successfully applied to solve these ambiguity problems. Semantic similarity quantifies how close two classes are in terms of semantics encoded in a given ontology. Using the web tool Semantic Similarity Measures using Disjunctive Shared Information ([DiShIn](http://labs.rd.ciencias.ulisboa.pt/dishin/)), we can calculate the semantic similarity between our recognized classes. For example, we can calculate the similarity between _Kaufman oculocerebrofacial syndrome_
(DOID:0111456) and _Syndromic intellectual disability_ (DOID:0050888), and the similarity between _Kagami-Ogata syndrome_ (DOID:0111712)
and _Syndromic intellectual disability_ (DOID:0050888).

We would see that for all measures _Syndromic intellectual disability_ is much more similar to _Kaufman oculocerebrofacial syndrome_ than to _Kagami-Ogata syndrome_. This means that by using semantic similarity we can automatically identify _Kaufman oculocerebrofacial syndrome_ as the correct linked entity for the mention KOS in this text.

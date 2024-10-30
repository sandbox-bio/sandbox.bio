<script>
import Execute from "$components/Execute.svelte";
</script>

## Relation Extraction

Finding the relevant entities in text is sometimes not enough. We need to
know which sentences may describe possible relationships between those entities, such as a relation between a disease and a compound.

This a complex text mining challenge, but a simple approach is to construct a pattern that allow any kind of characters between two entities:

<Execute command="grep -n -w -E 'MH[SNE]?.*(C|c)affeine' chebi_27732_sentences.txt" />

The following sentence is one of the eight displayed sentences mentioning
a possible relation:

```text
... MHS families were investigated 
    with a caffeine ...
```

However, we are missing all the sentences that have caffeine first:
<Execute command="grep -n -w -E '(C|c)affeine.*MH[SNE]?' chebi_27732_sentences.txt" />

We will be able to see that sometimes caffeine comes first:

```text
... caffeine-halothane contracture test were
    greater in those who had a known MH ...
... caffeine threshold and tension values
    and the MH ...
```

#### Multiple filters

The most flexible approach is use two `grep` commands. The first selects the
sentences mentioning one of the entities, and the other selects from the previously selected sentences the ones mentioning the other entity. For example,
we can first search for the acronyms and then for caffeine:

<Execute command="grep -n -w -E 'MH[SNE]?' chebi_27732_sentences.txt | grep -w -E '(C|c)affeine'" />

This will show all the ten sentences mentioning caffeine and an acronym.

#### Relation type

If we are interested in a specific type of relationship, we may have an additional filter for a specific verb. For example, we can add a filter for sentences with the verb response or diagnosed:

<Execute command="grep -n -w -E 'MH[SNE]?' chebi_27732_sentences.txt | grep -w -E '(C|c)affeine' | grep -w -E 'response|diagnosed'" />

We should note that this does not take in account where the verb appears
in the sentence. For example, in the following sentence the verb response
appears first than any of the two entities:

```text
The relationship between the IVCT response
and genotype was ... the number of MHS discordants
... at 2.0 mM caffeine ...
```

If the verb needs to appear between the two entities, we have to construct
a pattern that have these words in the middle of them:

<Execute command="grep -n -w -E 'MH[SNE]?.*(response|diagnosed).*(C|c)affeine' chebi_27732_sentences.txt" />

We can see now that the previous sentence is not presented as a match.

#### Remove relation types

We may also be interested in ignoring specific type of relations. To do that,
we only need to use the -v (or --invert-match) option. For example, to
ignore sentences with the word response or diagnosed:

<Execute command="grep -n -w -E 'MH[SNE]?' chebi_27732_sentences.txt | grep -w -E '(C|c)affeine' | grep -v -w -E 'response|diagnosed'" />

All the resulting sentences do not mention response or diagnosed.

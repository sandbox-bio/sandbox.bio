<script>
import Execute from "$components/Execute.svelte";
</script>

## Pattern File

The script created in the previous section only accepts one pattern, however
we may need to recognize different entities, or different mentions of the
same entity, such as the official name, possible synonyms, and the acronyms.
Fortunately, `grep` allows us to include a list of patterns directly from a file
using the `-f` option. The equivalent long form to the `-f` option is `--file=FILE`. 
For example, we can create a text file named `patterns.txt` 

<Execute command="nano patterns.txt" />

with the following three patterns:
```bash
(M|m)alignant (H|h)yperthermia
MH[SNE]?
(C|c)affeine
```

Then we can execute the previous grep but using multiple patterns specified in the pattern file:

<Execute command="grep -n -o -w -E -f patterns.txt chebi_27732_sentences.txt" />

Analyzing the output, we can check that the same sentences may include
different entities.
We can now update our script named `getentities.sh` 

<Execute command="nano getentities.sh" />

to receive as input not a single pattern but the filename where multiple patterns can be found.

```bash
PATTERNS=$1
grep -n -o -w -E -f $PATTERNS | \
tr ':' '\t'
```

We can execute the script giving as argument the file containing the patterns:

<Execute command="./getentities.sh patterns.txt < chebi_27732_sentences.txt" />

To save the output as a file named chebi_27732.tsv, we only need to add
the redirection operator: 

<Execute command="./getentities.sh patterns.txt < chebi_27732_sentences.txt > chebi_27732.tsv" />

Using the `patterns.txt` file is very useful if for example we are not focused
in a single disease, and we want to find any disease mentioned in the text.
In these cases, we have to create a file with the full lexicon of diseases. This
topic will be addressed in the following tutorial.
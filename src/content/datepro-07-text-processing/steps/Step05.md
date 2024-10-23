<script>
import Execute from "$components/Execute.svelte";
</script>

## Entity recognition

To select the sentences with one of our acronyms, we can use the grep command and our sentences file:

<Execute command="grep -w -E 'MH[SNE]?' chebi_27732_sentences.txt" />

The output will only include matching sentences. 

Alternatively, we can use the `-n` option to get the number of the line and
the `-o` option to get the acronym matched :

<Execute command="grep -n -o -w -E 'MH[SNE]?' chebi_27732_sentences.txt" />

The equivalent long form to the `-n` option is `--line-number`. 

We can also add the `-b` option to get the exact position of the acronym
matched:

<Execute command="grep -b -n -o -w -E 'MH[SNE]?' chebi_27732_sentences.txt" />

The output now contains the number of the line, the character position,
and the match.

We can now make a script that receives a pattern as argument and the
input text as the standard input, to display the line numbers and the matches
in a TSV format. Thus, let us create a script file named `getentities.sh`:

<Execute command="nano getentities.sh" />

with the following lines:
```bash
PATTERN=$1
grep -n -o -w -E $PATTERN | \
tr ':' '\t'
```
Again we should not forget to save the file in our working directory, and add
the right permissions with `chmod`, as we did with our scripts in the previous tutorials.

<Execute command="chmod u+x getentities.sh" />

The first line stores the pattern given as argument in the variable `PATTERN`. The grep command finds the matches and the `tr` command replaces each
colon by a tab character to produce TSV content.
We can now execute the script giving the pattern as argument and the
sentences file as standard input:

<Execute command="./getentities.sh 'MH[SNE]?' < chebi_27732_sentences.txt" />

We should note that now we have the values separated by a tab character,
i.e. the output is in TSV format.
The output can also be saved as a TSV file that we can open directly in our
preferred spreadsheet application. 
For example, to save it as chebi_27732.tsv, we only need to add the redirection operator:
<Execute command="./getentities.sh 'MH[SNE]?' < chebi_27732_sentences.txt > chebi_27732.tsv" />

#### Select the sentence
If we want to analyze a specific matched sentence, we can use a text editor
and go to that line number. A more efficient alternative is to use the print `p`
option of `sed` to output a given line number. For example, to check the MHS
match at line 3:

<Execute command="sed -n '3p' chebi_27732_sentences.txt" /> 

Now we can easily check the context of the match.
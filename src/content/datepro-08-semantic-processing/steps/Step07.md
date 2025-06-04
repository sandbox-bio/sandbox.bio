<script>
import Execute from "$components/Execute.svelte";
</script>

## Generic Lexicon

Instead of using a customized and limited lexicon, we may be interested in recognizing any of the diseases represented in the ontology. By recognizing all the diseases in our caffeine related text, we will be able to find all the diseases that may be related to caffeine

#### All labels

To extract all the labels from the disease ontology we can use the same XPath query used before, but now without restricting it to any URI:

<Execute command={`xmllint --xpath "//*[local-name()='Class']/*[local-name()='hasExactSynonym' or local-name()='hasRelatedSynonym' or local-name()='label']/text()" doid.owl`} />

We can create a script named `getalllabels.sh`:

<Execute command="nano getalllabels.sh" />

That receives as argument the OWL file where to find all labels containing the following lines:

```bash
OWLFILE=$1

xmllint --xpath "//*[local-name()='Class']/*[local-name()='hasExactSynonym' or local-name()='hasRelatedSynonym' or local-name()='label']/text()" $OWLFILE | sort -u
```

We should note that this script is similar to the `getlabels.sh` script without the `xargs`, since it does not receive a list of URIs as standard input.
Now we can execute the script to extract all labels from the OWL file:

<Execute command="chmod u+x getalllabels.sh" />

<Execute command="./getalllabels.sh doid.owl" />

The output will contain the full list of diseases in this reduced version of the ontology.

To create the generic lexicon, we can redirect the output to the file diseases.txt:

<Execute command="./getalllabels.sh doid.owl > diseases.txt" />

We can check how many labels we got by using the wc command:

<Execute command="wc -l diseases.txt" />

From the original OWL file the lexicon contains more than 34 thousand labels.

We can now recognize the lexicon entries in the sentences of the file
`chebi_27732_sentences.txt` by using the grep command:

<Execute command="grep -n -w -F -f diseases.txt chebi_27732_sentences.txt" />

The output will show the large list of sentences mentioning diseases.

Problematic entries
Despite using the `-F` option, the lexicon contains some problematic entries.
Some entries have expressions enclosed by parentheses or brackets, that represent alternatives or a category:

```text
Post measles encephalitis (disorder)
Glaucomatous atrophy [cupping] of optic disc
```

Other entries have separation characters, such as commas or colons, to
represent a specialization. For example:

```text
Tapeworm infection: intestinal taenia solum
Tapeworm infection: pork
Pemphigus, Benign Familial
ATR, nondeletion type
```

A problem is that not all have the same meaning. A comma may also be
part of the term. For example:

```text
46,XY DSD due to LHB deficiency
```

Other case includes using `&amp;` to represent an ampersand. For example:

```text
Gonococcal synovitis &amp;/or tenosynovitis
```

However, most of the times the alternatives are already included in the
lexicon in different lines. For example:

```text
Gonococcal synovitis and tenosynovitis
Gonococcal synovitis or tenosynovitis
```

As we can see by these examples, it is not trivial to devise rules that fully
solve these issues. Very likely there will be exceptions to any rule we devise and that we are not aware of.

#### Special characters frequency

To check the impact of each of these issues, we can count the number of times they appear in the lexicon:

<Execute command="grep -c -F '(' diseases.txt" />

<Execute command="grep -c -F ',' diseases.txt" />

<Execute command="grep -c -F '[' diseases.txt" />

<Execute command="grep -c -F ':' diseases.txt" />

<Execute command="grep -c -F '&amp;' diseases.txt" />

In the original OWL file we would see that parentheses and commas are the most frequent, with more than one thousand entries.

#### Completeness

Now let us check if the ATR acronym representing the _alpha thalassemia-X- linked intellectual disability syndrome_ is in the lexicon:

<Execute command="grep -E '^ATR' diseases.txt" />

All the entries include more terms than only the acronym.

Thus, a single ATR mention will not be recognized.
This is problematic if we need to match sentences mentioning that acronym, such as:

<Execute command="echo 'The ATR syndrome is an alpha thalassemia that has material basis in mutation in the ATRX gene on Xq21' | grep -w 'ATR'" />

We will now try to mitigate these issues as simply as we can. We will not
try to solve them completely, but at least address the most obvious cases.

#### Removing special characters

The first fix we will do, is to remove all the parentheses and brackets by using the tr command, since they will not be found in the text:

<Execute command="tr -d '[]()&lcub;&rcub;' < diseases.txt" />

Of course, we may lose the shorter labels, such as _Post measles encephalitis_, but at least now, the disease _Post measles encephalitis disorder_ will be recognized:

<Execute command="tr -d '[]()&lcub;&rcub;' < diseases.txt | grep 'Post measles encephalitis disorder'" />

If we really need these alternatives, we would have to create multiple
entries in the lexicon or transform the labels in regular expressions.

#### Removing extra terms

The second fix is to remove all the text after a separation character, by using the `sed` command:

<Execute command="tr -d '[]()&lcub;&rcub;' < diseases.txt | sed -E 's/[,:;] .*$//'" />

We should note that the regular expression enforces a space after the separation character to avoid separation characters that are not really separating two expressions, such as: _46,XY DSD due to LHB deficiency_
We can see that now we are able to recognize both ATR and _ATR syndrome_:

<Execute command="tr -d '[]()&lcub;&rcub;' < diseases.txt | sed -E 's/[,:;] .*$//' | grep -E '^ATR'" />

#### Removing extra spaces

The third fix is to remove any leading or trailing spaces of a label:

<Execute command="tr -d '[]()&lcub;&rcub;' < diseases.txt | sed -E 's/[,:;] .*$//; s/^ *//; s/ *$//'" />

We should note that we added two more replacement expressions to the `sed` command by separating them with a semicolon.

We can now update the script `getalllabels.sh`:

<Execute command="nano getalllabels.sh" />

To include the previous `tr` and `sed` commands:

```bash
OWLFILE=$1

xmllint \
    --xpath \
        "//*[local-name()='Class']/*[local-name()='hasExactSynonym' or local-name()='hasRelatedSynonym' or local-name()='label']/text()" $OWLFILE | \
    tr -d '[](){}' | \
    sed -E 's/[,:;] .*$//; s/^ *//; s/ *$//' | \
    sort -u
```

And we can now generate a fixed lexicon:

<Execute command="./getalllabels.sh doid.owl > diseases.txt" />

We can check again the number of entries:

<Execute command="wc -l diseases.txt" />

We have less entries because our fixes made some entries equal to others already in the lexicon, and thus the `-u` option filtered them.
From the original OWL we would get a lexicon with more than 13 thousand labels.

#### Disease recognition

We can now try to recognize lexicon entries in the sentences of file `chebi_27732_sentences.txt`:

<Execute command="grep -n -o -w -F -f diseases.txt chebi_27732_sentences.txt" />

To obtain the list of labels that were recognized, we can use the grep
command:

<Execute command="grep -o -w -F -f diseases.txt chebi_27732_sentences.txt | sort -u" />

From the original OWL we would get a list of 47 unique labels representing diseases that may be related to caffeine:

The reason why `47` appears is because there is a label 47, XXY:

<Execute command="echo '47, XXY' | ./geturi.sh doid.owl" />

#### Case insensitive

We may use the `-i` option to perform a case insensitive matching. To check how many labels are now being recognized we can execute:

<Execute command="grep -o -w -F -i -f diseases.txt chebi_27732_sentences.txt | sort -u | wc -l" />

From the original OWL we would get 66 labels being recognized.

To check which new labels were recognized, we can compare the results
with and without the `-i` option:

<Execute command="grep -o -w -F -i -f diseases.txt chebi_27732_sentences.txt | sort -u > diseases_recognized_ignorecase.txt" />

<Execute command="grep -o -w -F -f diseases.txt chebi_27732_sentences.txt | sort -u > diseases_recognized.txt" />

<Execute command="grep -v -F -f diseases_recognized.txt diseases_recognized_ignorecase.txt" />

We are now able to see that the new labels are.

Some of them are just lower and upper case variations of the same label.

To verify this, we can add the `-f` option to the sort command:

<Execute command="grep -o -w -F -i -f diseases.txt chebi_27732_sentences.txt | sort -u -f | wc -l" />

From the original OWL we would get 57 different labels being recognized. The equivalent long form to the `-f` option is `--ignore-case`.
Correct matches

Some important diseases could only be recognized by performing a case insensitive match, such as dyskinesia. This disease was missing because in the lexicon we had the uppercase case version of the labels, but not the lowercase version. We can check it by using the grep command:

<Execute command="grep -i -E '^dyskinesia$' diseases.txt" />

The lexicon has only the disease name with the first character in uppercase.

#### Incorrect matches

However, using a case insensitive match may also create other problems, such as the acronym CAN for the disease _Crouzon syndrome-acanthosis nigricans syndrome_:

<Execute command="echo 'CAN' | ./geturi.sh doid.owl | ./getlabels.sh doid.owl" />

By using a case insensitive grep we will recognize the common word CAN as a disease. For example, we can check how many times CAN is recognized:

<Execute command="grep -n -o -w -i -F -f diseases.txt chebi_27732_sentences.txt | grep -i ':CAN' | wc -l" />

And to see which type of matches they are, we can execute the following command:

<Execute command="grep -o -w -i -F -f diseases.txt chebi_27732_sentences.txt | grep -i -E '^CAN$' | sort -u" />

We can verify that the matches are incorrect mentions of the disease
acronym: `can` This means we created many mismatches by performing a case insensitive match.

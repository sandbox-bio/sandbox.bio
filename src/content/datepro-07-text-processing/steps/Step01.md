<script>
import Execute from "$components/Execute.svelte";
</script>

## Pattern Matching 

We used the `grep` command in the last tutorials to find a disease in the text,
since grep receives as argument a pattern to find an exact match in the text,
like any search functionality provided by conventional text editors. However,
we may need to search for multiple patterns even when interested in a single
disease. For example, when searching for mentions of malignant hyperthermia, we may also be interested in finding mentions using related expressions,
such as:
- MH : acronym
- MHS : acronym for malignant hyperthermia susceptible

Since we already know how to deal with multiple patterns by using the `-e` option, we may easily solve this problem by executing:

<Execute command="grep -e 'malignant hyperthermia' -e 'MH' -e 'MHS' chebi_27732.txt" />

#### Case insensitive matching

When dealing with text, using a case sensitive search is usually a good approach to avoid wrong matches. For example, acronyms are normally in upper case, while the full name is usually in lowercase having sometimes the
first letter of each word (or only the first word) in uppercase. So, instead of
using a full case sensitive grep, we might think on performing a case sensitive `grep` for the acronyms and a case insensitive grep for the disease words using the `-i` option:

<Execute command="grep -e 'MH' -e 'MHS' chebi_27732.txt" />

<Execute command="grep -i -e 'malignant hyperthermia' chebi_27732.txt"/>

The equivalent long form to the `-i` option is `--ignore-case`. We should note that each execution of grep will produce two separate lists of matching lines that might be overlapped. Alternatively, we can also convert it to just one case sensitive grep, if we are sure that Malignant hyperthermia is the only alternative case to malignant
hyperthermia present in the text. So, we can add it as another pattern:

<Execute command="grep -e 'Malignant hyperthermia' -e 'malignant hyperthermia' -e 'MH' -e 'MHS' chebi_27732.txt" /> 

#### Number of matches
To be sure that we are not losing any match, we can count the number of matching lines for both cases. First we execute a case insensitive grep and then we execute a case sensitive `grep`, both using the `-c` option:

<Execute command="grep -c -i 'malignant hyperthermia' chebi_27732.txt" />

<Execute command="grep -c -e 'malignant hyperthermia' -e 'Malignant hyperthermia' chebi_27732.txt" />

The equivalent long form to the `-c` option is `--count`.
In our case, the output should show 100 and 98 matching lines for the
insensitive and sensitive patterns, respectively.
This means that there is two lines that were not caught by the case sensitive pattern. To identify them, we can manually analyze each of the 100 matching lines one by one. But the goal of this tutorial is exactly avoiding these type of tedious tasks. One thing we can do to solve this issue is to find from the case insensitive matches the one that do not match the case sensitive patterns.

### Invert match

Fortunately, the grep command has the `-v` option that inverts the matching
and returns the lines of text that do not contain any matching. The equivalent
long form to the `-v` option is `--invert-match`.
Thus, if we apply the inverted match with the case sensitive patterns to
the output given by the case insensitive matching, we will get our outlier
mention:

<Execute command="grep -i 'malignant hyperthermia' chebi_27732.txt | grep -v -e 'Malignant hyperthermia' -e 'malignant hyperthermia'" />

From the output, we can easily identify the missing matching lines:

```text
...gene are associated with Malignant Hyperthermia (MH) and...
```

We were missing the case where both words have the first letter in uppercase.
Thus, to obtain all the matching lines in a case sensitive match we just
have to include the missing match as another pattern:

<Execute command="grep -c -e 'malignant hyperthermia' -e 'Malignant hyperthermia' -e 'Malignant Hyperthermia' chebi_27732.txt" />

#### File Differences

Another alternative to compare different matches, is to use the diff command that receives as input two files and identifies their differences. So, we
can create two auxiliary files and then apply the diff to them:

<Execute command="grep -i 'malignant hyperthermia' chebi_27732.txt > insensitive.txt" />

<Execute command="grep -e 'Malignant hyperthermia' -e 'malignant
hyperthermia' chebi_27732.txt > sensitive.txt" />

<Execute command="diff sensitive.txt insensitive.txt" />

The output should be the same text.
Alternatively, we can use process substitution that allows the output of a
command to be used as a file input to another command. This way, we can
apply the two grep commands directly in diff without needing extra files:

<Execute command="diff <(grep -i 'malignant hyperthermia' chebi_27732.txt) <(grep -e 'Malignant hyperthermia' -e 'malignant hyperthermia' chebi_27732.txt)" />

A problem that may occur with case sensitive matching is that some
acronyms are defined with lowercase letters in the middle, such as ChEBI,
and humans are not consistent with the way they mention them. The same
acronym may be mentioned in their original form or with all letters in upper-
case, or just some of them. Moreover, these inconsistent mentions sometimes
may even be found in the same publication.

#### Evaluation metrics
These inconsistencies made by humans when mentioning case sensitive expressions, is one of the reasons that most online search engines use case insensitive searches as default. This type of approach favors recall, while case sensitive search favor precision. Recall is the proportion of the number of correct matches found by our tool over the total number of correct mentions in the texts (found or not found).
Case insensitive searches avoid missing mentions, so they favor recall.
Precision is the proportion of the number of correct matches found by
our tool over the total number of matches found (correct or incorrect). Case
sensitive searches avoid incorrect matches, so they favor precision.
Normally, there is a trade-off between precision and recall. Using a tech-
nique that improves precision, most of the times, will decrease recall, and
vice-versa. To know how good the trade-off is, we can use the [F-measure](https://en.wikipedia.org/wiki/F1_score),
which is the harmonic average of the [precision and recall](https://en.wikipedia.org/wiki/Precision_and_recall).

#### Word Matching
Acronyms (or terms) may also appear inside common words or longer
acronyms. For example, when searching for MH, the word _victimhood_ will
produce a match:

<Execute command={`$echo "victimhood" | grep -i 'MH'`} />

The problem with _victimhood_ could be easily solved by using case sensitive
matching, but not for a longer acronym. For example, the acronym NEDMHM
for neurodevelopmental disorder with midbrain and hindbrain malformations
will produce a case sensitive match:

<Execute command={`echo "NEDMHM" | grep 'MH'`} />

One way to address this problem is to use the -w option of grep to only
match entire words, i.e. the match must be preceded and followed by char-
acters that are not letters, digits, or an underscore (or be at the beginning or end of the line). The equivalent long form to the `-w` option is `--word-regexp`.

Using this option, neither victimhood or NEDMHM will produce a match:

<Execute command={`echo "victimhood" | grep -w -i 'MH'`} />
<Execute command={`echo "NEDMHM" | grep -w -i 'MH'`} />

Word matching improves precision but decreases recall, since we may miss
some less common acronyms that we are not aware of, but are still relevant
for our study. For example, consider that we may also be interested in the
following acronyms:
- MHE : acronym for malignant hyperthermia equivocal
- MHN: acronym for malignant hyperthermia normal

If we apply word matching, we will not get a match, since both exact
matches are followed by a letter:

<Execute command={`echo "MHE and MHN" | grep -w -i 'MH'`} />

These are not trivial problems to solve by exact pattern matching, we may
need regular expressions to address some of these issues more efficiently.

<script>
import Execute from "$components/Execute.svelte";
</script>

## Position

Sometimes besides the match, we are also interested in limiting the matches
to specific parts of the input string. For example, to identify start and stop
codons in a protein sequence, we need to limit the matches to the beginning
or the end of the sequence. In text, we may for example be interested in
lines starting with a name of a disease. To take in account the position of a
match, regular expressions patterns can start with the circumflex character
(`^`) and/or end with the dollar sign character (`$`).
If the pattern starts with a circumflex then only matches at the beginning
of the line will be considered. On the other hand, if the pattern ends with a
dollar then only matches at the end of the line will be considered.

#### Beginning

For example, if we are looking for lines starting with Malignant Hyperthermia
we can use the following pattern:

<Execute command="grep -E '^(M|m)alignant (H|h)yperthermia' chebi_27732.txt" />

The output will include the list of lines beginning with a mention to Malignant Hyperthermia:

```text
... Malignant hyperthermia (MH) is a
potentially fatal autosomal ...

Malignant hyperthermia (MH) is a
pharmacogenetic disorder ...
```

To check how many of the matching lines were filtered, we can count the
number of occurrences when using the circumflex and when not:

<Execute command="grep -c -E '^(M|m)alignant (H|h)yperthermia' chebi_27732.txt
grep -c -E '(M|m)alignant (H|h)yperthermia' chebi_27732.txt" />

The output will show that only 20 of the 100 matches were considered.

#### Ending

If we are looking for lines ending with a mention to Malignant Hyperthermia,
then we can add the dollar character to the end of the pattern:

<Execute command="grep -E '(M|m)alignant (H|h)yperthermia.$' chebi_27732.txt" />

To allow a punctuation character before the end of the line, we added the
dot character before the dollar character in the pattern. The dot character
matches any character, including the dot itself.

The output will be the list of lines ending with a mention to Malignant
Hyperthermia.

We can check how many lines were filtered by using again the -c option:
<Execute command="grep -c -E '(M|m)alignant (H|h)yperthermia.$' chebi_27732.txt" />
<Execute command="grep -c -E '(M|m)alignant (H|h)yperthermia' chebi_27732.txt" />

The output will show that only 15 of the 100 matches were at the end of the
line.

#### Near the end

Sometimes we do not want the mention ending exactly at the last character.
We may be more flexible and allow a following expression, or a given number
of characters. For example, to allow 10 other characters between the end of
the line and the mention of Malignant Hyperthermia, we can add a quantifier
to the dot operator:

<Execute command="grep -c -E '(M|m)alignant (H|h)yperthermia.{0,10}$' chebi_27732.txt" />

The output will show that we have 20 matches.
If we remove the `-c` option, we will be able to check that words, such
as families and patients, are now allowed to appear between the mention of
Malignant Hyperthermia and the end of the line.

#### Word in between

To allow a word in between, independently of its length, we can add to the
pattern an optional sequence of non-space characters (the word) preceded
by a space:

<Execute command="grep -c -E '(M|m)alignant (H|h)yperthermia( [^ ]*)?.$' chebi_27732.txt" />

The output will show that we have 24 matches. We should note that the `[^ ]`
operator avoids having two words.

If we remove the `-c` option, we will be able to check that lengthy words
(with more than 10 characters), such as _susceptibility_, are now allowed to
appear between the mention of Malignant Hyperthermia and the end of the
line.

##### Full line

If we want lines that start with a mention to Malignant Hyperthermia and end
with an acronym, MH or MHS, then we can execute two grep commands.
The first gets the lines starting with Malignant Hyperthermia and the next
filters the output of the latter with lines ending with an acronym:

<Execute command="grep -E '^(M|m)alignant (H|h)yperthermia' chebi_27732.txt | grep -w -E 'MHS?.$'" />

Alternatively, we can add both the circumflex and dollar operators to the
same pattern. However, we cannot forget to add `.*` to match anything in
between them, since we are asking full line matches:

<Execute command="grep -w -E '^(M|m)alignant (H|h)yperthermia.*MHS?.$' chebi_27732.txt" />

We can see that both commands match all the text of the abstract since
each abstract is stored in a single line of the file.

This demonstrates the problem of tokenization, since usually what we really
need is to match a full sentence or a phrase. And in that case each line should
represent a sentence or phrase from the abstract.

#### Match position

For more advanced processing, we may be interested in knowing the exact
position of the matches in a given line. This can be done by using the `-b`
option of grep, which provides the number of bytes in the line before the
start of the match:

<Execute command="echo 'MHS MHN MHE' | grep -b -o -w -E 'MH[SNE]'" />

The equivalent long form to the `-b` option is `--byte-offset`.
The output shows the list of matches preceded by their position.

The same result happens if the input is given in multiple lines:

<Execute command="echo -e 'MHS\nMHN\nMHE' | grep -b -o -w -E 'MH[SNE]'" />

We have the exact same result because the newline character counts the same
as the space.

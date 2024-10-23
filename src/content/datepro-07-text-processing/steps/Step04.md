<script>
import Execute from "$components/Execute.svelte";
</script>

## Tokenization

As we have shown in the previous section, sometimes we need to work at the
level of a sentence and not use a full document as the input string. Tokeniza-
tion is a Natural Language Processing (NLP) task that aims at identifying
boundaries in the text to fragment it into basic units called tokens. These
tokens can be sentences, phrases, multi-word expressions, or words.

#### Character delimiters

In most languages, some specific characters can be considered as accurate
boundaries to fragment text into tokens. For example, the space character
to identify words; the period (.), the question mark (?) and the exclamation
mark (!) to identify the ending of a sentence; and the comma (,), the semi-
colon (;), the colon (:) or any kind of parenthesis to identify a phrase within a
sentence. However, this problem may be more complex in languages without
explicitly delimiters, such as Chinese.

A common approach to tokenization is to use regular expressions to replace these delimiters by newline characters. This will result in a token per-line. For example, we can replace the characters specifying the end of a sentence with a newline by using the tr command and then count the number of lines:

<Execute command="tr '[.!?]' '\n' < chebi_27732.txt | wc -l" />

We get 1618 lines from the original 255 lines.

<Execute command="wc -l chebi_27732.txt" />

Unfortunately, this is not just so simple. We need to analyze the output:

<Execute command="tr '[.!?]' '\n' < chebi_27732.txt | less" />

#### Wrong tokens

We can check that: i) many lines are empty because an extra newline character will be added to the last sentence, and ii) the dot character is also used
as a decimal mark in a number, then some sentences are split in multiple
lines because they have decimal number in them. For example, the original
sentence:
```text
These 10 mutations account for 21.9% of the North American MH-susceptible population
```
is split in two lines:
```text
These 10 mutations account for 21
9% of the North American MH-susceptible population
```

#### String Replacement
This means that looking at just one character is not enough, we need some
context. For performing this, we will use the `sed` command that we may
consider as a more powerful version of the `tr` command. The `sed` command
is a stream editor that can receive as input a string and perform basic text
transformations, such as replace one expression by another, that are available
in almost all text editors. For example, we can use a simple `sed` to convert
every mention of caffeine by its ChEBI identifier:

<Execute command="sed -E 's/caffeine/CHEBI:27732/gi' chebi_27732.txt" />

The `-E` option allow us to use extended regular expressions, like we used
before in `grep`. The `s` option has the following syntax `'s/FIND/REPLACE/FLAGS'`, where: `FIND` is the pattern to find in the input string; `REPLACE`
the expression to replace the matches; `FLAGS` are multiple options, such as
`g` to replace all matches in each line and not just the first one, and `i` to be
case insensitive.
For example, the original fragment of text:
``text
... link between the caffeine threshold and tension ...
```
will be converted to:
```text
... link between the CHEBI:27732 threshold and tension ...
```

#### Multi-character delimiters
To replace the delimiter characters by a newline when followed by at least
one space character, we can use the following command:

<Execute command="sed -E 's/[.!?] +/\n/g' chebi_27732.txt" />

We should note that by making compulsory a space character, we avoid: i)
empty lines by splitting a sentence that is already at the end of the line (assuming there are no ghost space characters at the end of each line), and ii) decimal markers because they are followed by numerical digits and not spaces.
We now get 1092 lines from the original 255 lines:

<Execute command="sed -E 's/[.!?] +/\n/g' chebi_27732.txt | wc -l" />

#### Keep delimiters
The previous `sed` command is removing the delimiter characters from the
text, and this may cause other problems. A better solution is to keep the
delimiter characters and just add the newline. The `sed` command allows us to
keep each match for a specific part of the pattern (sub-pattern) by enclosing
it within parentheses. To include the match of a sub-pattern in the replace
expression, we can use the backslash and its numerical order. Thus, we can
improve our `sed` command by using this technique so we do not remove any
delimiter character:

<Execute command="sed -E 's/([.!?])( +)/\1\n\2/g' chebi_27732.txt" />

The `\1` represents the match for the sub-pattern `([.!?])`, and the `\2` represents the match for the sub-pattern (` +`). This means that a newline character is inserted right after each delimiter character found, and keeping the space characters.

For example, the original fragment of text:

```text
... muscle relaxants. To date, ...
```
will be converted to:
```text
... muscle relaxants.
To date, ...
```

However, other common issues may still persist:
```text
... bulk.&lt;h4&gt;Methods&lt;/h4&gt;Fetal ...
... sequencing.&lt;h4&gt;Results&lt;/h4&gt;Whole ...
```

These sentences include a HTML elements.
To minimize this issue, we can change the pattern to add the option of `&`
character besides the space:

<Execute command="sed -E 's/([.!?])([& ]+)/\1\n\2/g' chebi_27732.txt | wc -l" />

We now get 1179 lines, i.e. this pattern is more flexible and was able to split
more 87 sentences. 

<Execute command="expr 1179 - 1092" />

This does not mean that is free of errors. It is almost impossible to derive
a rule that covers all the possible typos humans can produce.

#### Sentences file

To save the output as a file named `chebi_27732_sentences.txt`, we only need
to add the redirection operator:

<Execute command="sed -E 's/([.!?])([& ]+)/\1\n\2/g' chebi_27732.txt > chebi_27732_sentences.txt" />

Each line of the file `chebi_27732_sentences.txt` represents a sentence.
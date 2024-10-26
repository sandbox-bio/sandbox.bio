<script>
import Execute from "$components/Execute.svelte";
</script>

## Regular Expressions

When dealing with natural language text we may need more flexibility than
the one provided by exact matching. Regular expressions are an efficient
tool to extend exact matching with flexible patterns, that may find different
matches. As an example, we may be interested in finding all the mentions
of the acronym MHS or MHN in a text. For doing that, regular expressions
provide the alternation operator that helps us to solve this issue easily by
specifying multiple alternatives to match in a specific part of the pattern, in
this case an S or an N as the last character.
Regular expressions can be better understood by clearly separating three
distinct components:

- input : any string where we want to find something
- pattern : a string that specifies what we are looking for
- match : a fragment of the input (a substring) where the pattern can be found

In our examples, the input is the text file `chebi_27732.txt`, but it can be the
amino acid sequences that we previously extracted from the UniProt file entries. Until now the pattern has represented an exact string to look for, where
each match is an exact replica of the pattern occurring at a given position of
the input string. When using regular expressions, the pattern contains special characters, whose purpose are not to directly match with the input but instead have a special meaning. These special characters represent operators
that specify which different types of strings we want to find in the input. For
example, strings that start with MH and end with S or an N. By using regular
expressions, the matches are not replicas of the pattern, they can be different
strings as long as they satisfy the specified pattern.

###### Extended syntax

The grep command allows us the possibility to include regular expression
operators in the input pattern. grep understands two different versions of
regular expression syntax: [basic and extended](https://www.regular-expressions.info/posix.html). We will use the extended
syntax for two reasons: (i) the basic does not support relevant operators,
such as alternation; (ii) and to clearly differentiate exact matching from reg-
ular expression matching. Thus, we will start to use the `-E` option, which
makes the command interpret the pattern as an extended regular expression. The equivalent long form to the `-E` option is `--extended-regexp`.
We should note that this option does not affects the matching when using
a pattern without any regular expression operator, such as MH. For example,
the following commands will produce the same results:

<Execute command="echo -e 'MHS\nMHN' | grep 'MH'" />
<Execute command="echo -e 'MHS\nMHN' | grep -E 'MH'" />

Note, that we use the `-e` option so the echo command interpret the `\n`
characters as a newline. Thus, the echo command outputs two lines, that
are given as input to the grep command. We should note that the grep
command filters lines.

#### Alternation

The first regular expression operator we will test is the alternation, which
we introduced above. An alternation is represented by the bar character (`|`)
that specifies a pattern where any match must include either the preceding or
following characters. The preceding and following characters can be enclosed
within parentheses to better specify the scope of the alternation operator. For
example, the pattern for finding strings that start with MH and end with S or
an N can be written as:

<Execute command="echo -e 'MHS\nMHN' | grep -E 'MH(S|N)'" />

We can also use multiple patterns using the -E option:

<Execute command="echo -e 'MHS\nMHN' | grep -E -e 'MH(S|X)' -e 'MH(X|N)'" />

###### Basic syntax

If we use the basic regular expression syntax no match will be found, since
the alternation operator is not supported:

<Execute command=" echo -e 'MHS\nMHN' | grep 'MH(S|N)'" />

We will have a match only if the | and the parentheses are in the input
string, since it is not interpreted as an operator:

<Execute command="echo -e 'MH(S|N)' | grep 'MH(S|N)'" />

###### Scope

To better understand the scope of an alternation, we can remove the parentheses from the pattern and add the `-w` option:

<Execute command="echo -e 'MHS\nMHN' | grep -w -E 'MHS|N'" />

We only get the first line. This is explained because the alternation operator
is applied to all the preceding characters, i.e. the grep will search for the
MHS word or the N word. If we add a single N to the input string we already
get another match:

<Execute command="echo -e 'MHS\nN' | grep -w -E 'MHS|N'" />

We can also move the opening parenthesis one character to the left:

<Execute command="echo -e 'MHS\nMHN' | grep -E 'M(HS|N)'" />

Only MHS is now displayed, since the alternative now represents MN without
the H.

###### Multiple alternatives

We are not limited to two alternatives, we can have multiple `|` operators in
a pattern. For example, the following command will find any of the three
acronyms MHS, MHE or MHN:

<Execute command="echo -e 'MHS\nMHN\nMHE' | grep -E 'MH(S|N|E)'" />

We can now transform our previous grep command with multiple case
sensitive patterns:

<Execute command="grep -c -e 'Malignant hyperthermia' -e 'Malignant Hyperthermia' -e 'malignant hyperthermia' chebi_27732.txt" />

in a `grep` command with a single pattern using alternation:

<Execute command="grep -c -E '(M|m)alignant (H|h)yperthermia' chebi_27732.txt" />

And we will obtain the same 100 matching lines.

#### Multiple characters

A useful regular expression feature is that we can use the dot character (`.`) to
represent any character, so if we want to find all the acronyms that start with
MH we can execute the following command:

<Execute command="grep -o -w -E 'MH.' chebi_27732.txt | sort -u" />

We should note that we use the `-o` option of the command grep so it just displays the matches and not all the line that includes the match. The equivalent long form to the `-o` option is `--only-matching`.

The output will be the following three-character lines:

```text
MH
MH)
MH,
MH.
MH1
MH2
MHE
MHN
MHS
```

The `-o` option also solves the problem of counting the total number of
matches, and not just the number of lines with a match:

<Execute command="grep -o -w -E 'MH.' chebi_27732.txt | wc -l" />
<Execute command="grep -c -w -E 'MH.' chebi_27732.txt" />

The output will show that 164 matches were found in 47 lines. The `-c` option
overrides the `-o` option, i.e. if we use both in the same grep the output will
be just the number of lines.
If we really want to match only the dot character, we have to precede it
with a backslash character (`\`):

<Execute command="grep -o -w -E 'MH\.' chebi_27732.txt | sort -u" />

Now only the MH. will be displayed.
We can check that there are some matches that are not really acronyms,
such as MH) and MH,.

###### Spaces

We should note that MH appears because the space character can also be
matched. For example, the following text includes a word match with MH
since the parenthesis is considered a word delimiter character (not a letter, digit or underscore) :

```text
... susceptible to MH (MHS) ...
```

On the other hand, the following text does not include a word match with
MH :

```text
... markers and MH susceptibility ...
```

Thus, what we really want is matches where the third character is a letter or
a numerical digit.

Sometimes, the text includes other characters that also represent horizon-
tal or vertical space in typography, such as the tab character. All these characters are known as whitespaces and can be represented by the expression
`\s` in a pattern 4 . The following command demonstrates that both the space
and the tab characters are matched by `\s`:

<Execute command="echo -e 'space: :\ntab:\t:' | grep -E '\s'" />

###### Groups

Fortunately, the regular expressions include the group operator that let us
easily specify a set of characters. A group operator is represented by a set of
characters enclosed within square brackets. Any of the enclosed characters
can be matched.
For example, the previous command to find any of the three acronyms can
be replaced by:

<Execute command="echo -e 'MHS\nMHN\nMHE' | grep -E 'MH[SNE]'" />

We should note that only one of the three letters, S, N or E will be matched
in the input string.

###### Ranges

Still, this is not solving our need to only match letters or digits. However, we
can also specify characters ranges with the dash character (`-`). For example,
to find all the acronyms that start with MH followed by any alphabet letter:

<Execute command="grep -o -w -E 'MH[A-Z]' chebi_27732.txt | sort -u" />

This will result in only three acronyms.

We should note that `A-Z` represents any alphabet letter in uppercase, a
lowercase letter will not be matched:

<Execute command="echo -e 'MHS\nMHs' | grep -E 'MH[A-Z]'" />

The output will be only one line.

If we intend to keep the usage of a case sensitive grep and at the same time find lowercase matches, then we need to add the a-z range:

<Execute command="echo -e 'MHS\nMHs' | grep -E 'MH[A-Za-z]'" />

The output will be both lines.
We should note that the dot character inside a range represents itself and
not any character:

<Execute command="echo -e 'MHS\nMH.' | grep -E 'MH[.]'" />

The output will be only the last line.

Additionally, to include the acronyms that end with a numerical digit we
need to add the `0-9` range:

<Execute command="grep -o -w -E 'MH[A-Z0-9]' chebi_27732.txt | sort -u" />

Finally, we have the correct list of all three character acronyms starting
with MH.

###### Negation

Another frequent case is the need to match any character with a few exceptions. For example, if we need to find all the matches that start with MH followed by any character except an alphabet letter. Fortunately, we can use the negation feature within a group operator. The negation feature is represented
by the circumflex character (`^`) right next to the left bracket. The negation
means that all the characters and ranges enclosed within the brackets are the
ones that cannot be matched. Thus, a solution to the above example is to add
the `A-Z` range after the circumflex:

<Execute command="grep -o -w -E 'MH[^A-Z]' chebi_27732.txt | sort -u" />

We can see that all of the three acronyms MHS, MHE or MHN will be
missing from the output.
If we do not want the MH acronym, we can add the space character to
the negative group:

<Execute command="grep -o -w -E 'MH[^A-Z ]' chebi_27732.txt | sort -u" />

The output should now contain one less acronym.

#### Quantifiers

Above we were interested in finding acronyms composed of exactly three
characters. However, we may need to find all acronyms that start with MH
independently of their length. This functionality is also available in regular
expressions using the quantifiers operators.

###### Optional

The simplest quantifier is the optional operator that is specified by an item
followed by the question mark character (`?`). The item can be a character,
an operator or a sub-pattern enclosed by parentheses. That item becomes
optional for matching, i.e. a match can either contain that item or not.
For example, to find all the acronyms starting with MH and followed by
one alphabetic letter or none:

<Execute command="grep -o -w -E 'MH[A-Z0-9]?' chebi_27732.txt | sort -u" />

Given that the third character is optional the output will include the two-character acronym MH, but not the MH match.

We can add the space character to the group:

<Execute command="grep -o -w -E 'MH[A-Z0-9 ]?' chebi_27732.txt | sort -u" />

Now the output includes the two-character acronym MH and the MH match.

###### Multiple and optional

To find all the acronyms independently of their length, we can use the asterisk
character (`*`). The preceding item becomes optional and can be repeated
multiple times. For example, to find all the acronyms starting with MH and
which may be followed any number of alphabetic letters or numeric digits:

<Execute command="grep -o -w -E 'MH[A-Z0-9]*' chebi_27732.txt | sort -u" />

The output now includes the four-character acronym MHS1.

We should note that the grep command uses a greedy approach, i.e. it
will try to match as many characters as possible. For example, the following
command will match MH1 and not MH:

<Execute command="echo 'MH1' | grep -o -E 'MH[0-9]*'" />

####### Multiple and compulsory

To make the preceding item compulsory and able to repeat it multiple times,
we may replace the asterisk by the plus character (`+`). For example, the fol-
lowing pattern will find all the acronyms starting with MH followed by at
least one alphabetic letter or numeric digit:

<Execute command="grep -o -w -E 'MH[A-Z0-9]+' chebi_27732.txt | sort -u" />

We should note that the output does not contain the two character
acronym MH.

###### All options

The above quantifiers are the most popular, but the functionality of all of
them can be reproduced by using curly braces to specify the minimal and
maximum number of occurrences. The item is followed by an expression of
the type `&lcub;n,m&rcub;` where n and m are to be replaced by a number specifying the
minimum and maximum number of occurrences, respectively. n and m may
also be omitted, which means that no minimum or maximum limit is to be
imposed.

Using curly brackets, the question mark character (`?`) can be replaced by
`{0,1}`. Thus, the following two patterns are equivalent:

<Execute command="grep -o -w -E 'MH[A-Z0-9]?' chebi_27732.txt | sort -u" />

<Execute command="grep -o -w -E 'MH[A-Z0-9]{0,1}' chebi_27732.txt | sort -u" />

The asterisk character (`*`) can be replaced by `&lcub;0,&rcub;`. Thus, the following
two patterns are equivalent:

<Execute command="grep -o -w -E 'MH[A-Z0-9]*' chebi_27732.txt | sort -u" />

<Execute command="grep -o -w -E 'MH[A-Z0-9]&lcub;0,&rcub;' chebi_27732.txt | sort -u" />

The plus character (`+`) can be replaced by `&lcub;1,&rcub;`. Thus, the following two
patterns are equivalent:

<Execute command="grep -o -w -E 'MH[A-Z0-9]+' chebi_27732.txt | sort -u" />

<Execute command="grep -o -w -E 'MH[A-Z0-9]&lcub;1,&rcub;' chebi_27732.txt | sort -u" />

On the other hand using `{1,1}` is the same as not having any operator.
Thus, the following two patterns are equivalent:

<Execute command="grep -o -w -E 'MH[A-Z0-9]' chebi_27732.txt | sort -u" />

<Execute command="grep -o -w -E 'MH[A-Z0-9]{1,1}' chebi_27732.txt | sort -u" />

The previous commands display the all the three-character acronyms.
For example, if we are looking for acronyms with exactly 4 characters then
we can apply the following pattern:

<Execute command="grep -o -w -E 'MH[A-Z0-9]{2,2}' chebi_27732.txt | sort -u" />

We should note that we use 2 as both the minimum and maximum since MH
already count as 2 characters. The output of the previous command is now the four-character acronym: MHS1

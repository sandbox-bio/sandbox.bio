<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

#### Data Extraction

Some data in the CSV file may not be relevant regarding our information
need, i.e. we may need to identify and extract relevant data. In our case, we
will select the relevant proteins (lines) using the command line tool `grep`,
and secondly, we will select the column we need using the command line tool
`cut`. Since our information need is about diseases related to caffeine, we may
assume that we are only interested in proteins that have one of these topics
in the third column:

```text
CC - MISCELLANEOUS
CC - DISRUPTION PHENOTYPE
CC - DISEASE
```

Extracting lines from a text file is the main function of `grep`. The selection
is performed by giving as input a pattern that grep tries to find in each line,
presenting only the ones where it was able to find a match. The pattern is
the same as the one we normally use when searching for a word in our text
editor. The grep command also works with more complex patterns such as
regular expressions, that we will describe later on.

#### Single and multiple patterns

We can execute the following command that selects the proteins with the
topic `CC - MISCELLANEOUS`, our pattern, in our CSV file:

<Execute command="grep 'CC - MISCELLANEOUS' data/chebi_27732_xrefs_UniProt.csv" />

The `data` folder contains the files retrieved in the previous tutorial.

The output will be a shorter list of proteins, all with `CC - MISCELLANEOUS`
as topic.

The next step will show how to filter by multiple topics.

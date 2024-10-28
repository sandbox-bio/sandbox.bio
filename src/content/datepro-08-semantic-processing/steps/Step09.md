<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

## Large lexicons

The online tool [MER](https://github.com/lasigeBioTM/MER) is based on a shell script, so it can be easily executed as a command line to efficiently recognize and link entities using large lexicons.

#### MER installation

We preloaded a reduced version of MER with just the essential files. We can start by extracting it:

<Execute command="tar -xzf MER.tgz" />

Alternatively, if we want to run MER locally, we can download the latest compressed file (zip) version and extract its contents:

```bash
curl -O -L https://github.com/lasigeBioTM/MER/archive/
master.zip
unzip master.zip
mv MER-master MER
```

We now have to copy the Human Disease Ontology in to the data folder of MER, and then enter into the MER folder:

<Execute command="cp doid.owl MER/data/" />
<Execute command="cd MER" />

Lexicon files

To execute MER, we need first to create the lexicon files:

<Execute command="(cd data; ../produce_data_files.sh doid.owl)" />

This may take a few minutes to run. However, we only need to execute it
once, each time we want to use a new version of the ontology. If we wait, the output will include the last patterns of each of the lexicon files.
We can check the contents of the created lexicons by using the `tail` command:

<Execute command="tail data/doid_*" />

These patterns are created according to the number of words of each term.

#### MER execution

Now we are ready to execute MER, by providing each sentence from the file `chebi_27732_sentences.txt` as argument to its `get_entities.sh` script.

<Execute command={`cat ../chebi_27732_sentences.txt | tr -d "'" | xargs -I &lcub;&rcub; ./get_entities.sh '&lcub;&rcub;' doid`} />

We removed single quotes from the text, since they are special characters to the command line `xargs`. We should note that this is the get_entities.sh script inside the MER folder, not the one we created before. Now we will be able to obtain a large number of matches.

The first two numbers represent the start and end position of the match in the sentence. They are followed by the name of the disease and its URI in the ontology.

We can also redirect the output to a TSV file named `diseases_recognized.tsv`:

<Execute command="cat ../chebi_27732_sentences.txt | tr -d "'" | xargs -I &lcub;&rcub; ./get_entities.sh '&lcub;&rcub;' doid > ../diseases_recognized.tsv" />

We can now open the file in our spreadsheet application, such as LibreOffice Calc or Microsoft Excel.

<Alert>
By using the reduced OWL files preloaded in this tutorial, we will get just a portion of the full list of mentions of diseases from the Human Disease Ontology.
</Alert>

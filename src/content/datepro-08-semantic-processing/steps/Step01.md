<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

## Classes

In the previous tutorials we searched for mentions of caffeine and malignant
hyperthermia in text. However, we may miss related entities that may also
be of our interest. These related entities can be found in semantic resources,
such as ontologies. The semantics of caffeine and malignant hyperthermia are
represented in ChEBI and DO ontologies, respectively.

#### OWL files

The ontologies are preloaded for this tutorial, but not in their entirety:

<Execute command="du -h *.owl" />

As you can see, the preloaded files are only a few KB in size, while the original ones are several MB. Since the original files are very large OWL files containing numerous classes, we have reduced the number of classes to avoid long waiting times for command line execution.

<Alert>
Using the reduced OWL files preloaded in this tutorial, the command lines described here may work in a similar, but not identical, manner.
</Alert>

The original OWL files can be retrieved by using `curl`:

<pre class="code border p-2" style="white-space: pre-wrap">
curl -L -O "http://purl.obolibrary.org/obo/doid/releases/2021-03-29/doid.owl"

curl -L -O "http://purl.obolibrary.org/obo/chebi/198/chebi_lite.owl"</pre>

<Alert>The `curl` command here at this platform cannot download the files due to access restrictions. However, the same command will work just fine in your local terminal.</Alert>

The `-O` option saves the content to a local file named according to the name
of the remote file, usually the last part of the URL. The equivalent long form
to the `-O` option is `--remote-name`. The option `-L` enables the curl command to follow a [URL redirection](https://en.wikipedia.org/wiki/URL_redirection). The equivalent long form to the `-L`
option is `--location`.

The previous `curl` commands will create the files `chebi_lite.owl` and `doid.owl`, respectively.
We should note that these links are for the specific releases used in the book, and using another release may change the output of the examples presented in it.
To retrieve the most recent release we should use the following links:

- http://purl.obolibrary.org/obo/doid.owl
- http://purl.obolibrary.org/obo/chebi/chebi_lite.owl

To find other ontology links search for them on the [BioPortal](http://bioportal.bioontology.org/) or on the [OBO Foundry](http://www.obofoundry.org/) webpages. Alternatively, we can also get the OWL files from the [book file archive](http://labs.rd.ciencias.ulisboa.pt/book/).

#### Class label

Both OWL files use the XML format syntax. Thus, to check if our entities are
represented in the ontology, we can search for ontology elements that contain
them using a simple grep command:

<Execute command="grep '>malignant hyperthermia<' doid.owl" />

<Execute command="grep '>caffeine<' chebi_lite.owl" />

For each `grep` the output will be the line that describes the property label
(`rdfs:label`), which is inside the definition of the class that represents the
entity.

#### Class definition

To retrieve the full class definition, a more efficient approach is to use the
`xmllint` command, which we already used in previous tutorials to process XML:

<Execute command={`xmllint --xpath "//*[local-name()='label' and text()='malignant hyperthermia']/.." doid.owl`} />

The XPath query starts by finding the label that contains malignant hyperthermia and then `..` gives the parent element, in this case the Class element. From the output we can see that the semantics of malignant hyperthermia is much more than its label.

For example, we can check that malignant hyperthermia is a subclass of
(specialization) the entries [0050736](http://purl.obolibrary.org/obo/DOID_0050735) and [66](http://purl.obolibrary.org/obo/DOID_66). By clicking in the previous links will see
that malignant hyperthermia is a special case of a autosomal dominant disease
and of a muscle tissue disease.

We can search for those specific relations between malignant hyperthermia and the entries `0050736` and `66`:

<Execute command={`xmllint --xpath "//*[local-name()='label' and text()='malignant hyperthermia']/..//*[@*[local-name()='resource' and .='http://purl.obolibrary.org/obo/DOID_66' or .='http://purl.obolibrary.org/obo/DOID_0050736']]" doid.owl`} />

We added the `@*[local-name()='resource']` to extract the URI specified in an attribute resource of any descendant element `//*[...]`.
The relation specification uses the `subClassOf` element.

We can do the same to retrieve the full class definition of caffeine:

<Execute command={`xmllint --xpath "//*[local-name()='label' and text()='caffeine']/.." chebi_lite.owl`} />

From the output we can see that the types of semantics available for caffeine differs from the semantics of malignant hyperthermia, but they still share
many important properties, such as the definition of `subClassOf`.

The class caffeine is a specialization of two other entries: [26385](http://purl.obolibrary.org/obo/CHEBI_26385) (purine
alkaloid), and [27134](http://purl.obolibrary.org/obo/CHEBI_26385) (trimethylxanthine).
We can search for those specific relations between caffeine and the entries
`26385` and `27134`:

<Execute command={`xmllint --xpath "//*[local-name()='label' and text()='caffeine']/..//*[@*[local-name()='resource' and .='http://purl.obolibrary.org/obo/CHEBI_26385' or .='http://purl.obolibrary.org/obo/CHEBI_27134']]" chebi_lite.owl`} />

The relation specification uses the `subClassOf` element.

#### Related Classes

There are additional subclass relationships that do not represent subsump-
tion (_is-a_).
For example, the relationship between caffeine and the entry [25435](http://purl.obolibrary.org/obo/CHEBI_25435) (mutagen) is defined by the entry [0000087](http://purl.obolibrary.org/obo/RO_0000087) (has role) of the Relations Ontology.
This means that the relationship defines that caffeine has role mutagen.
We can search that specific relation between caffeine and mutagen (CHEBI:25435):

<Execute command={`xmllint --xpath "//*[local-name()='label' and text()='caffeine']/..//*[@*[local-name()='resource' and .='http://purl.obolibrary.org/obo/CHEBI_25435']]/../.." chebi_lite.owl`} />

The specification uses the Restriction element.

We can now search in the OWL file for the definition of the type of relation
has role (`RO:0000087`):

<Execute command={`xmllint --xpath "//*[local-name()='ObjectProperty'][@*[local-name()='about']='http://purl.obolibrary.org/obo/RO_0000087']" chebi_lite.owl`} />

The XPath query starts by finding the elements ObjectProperty and then
selects the ones containing the about attribute with the relation URI as
value.
We can check that the relation is neither transitive or cyclic.

<script>
import Execute from "$components/Execute.svelte";
</script>

## Parent Classes

Parent classes represent generalizations that may also be relevant to recognize in text. To extract all the parent classes of malignant hyperthermia, we can use the following XPath query:

<Execute command={`xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/DOID_8545']/*[local-name()='subClassOf']/@*[local-name()='resource']" doid.owl`} />

The first part of the XPath is the same as the above to get the class element, then `[local-name()='subClassOf']` is used to get the subclass
element, and finally `@*[local-name()='resource']` is used to get the attribute containing its URI.

The output should be the URIs representing the parents of class `8545`.

We can also execute the same command for caffeine:

<Execute command={`xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/CHEBI_27732']/*[local-name()='subClassOf']/@*[local-name()='resource']" chebi_lite.owl`} />

The output will now include the parents of caffeine.

We should note that we no longer can use the string function, because ontologies are organized as DAGs using multiple inheritance, i.e. each class can have multiple parents, and the string function only returns the first
match. To get only the URIs, we can apply the previous technique of using the `cut` command:

<Execute command={`xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/DOID_8545']/*[local-name()='subClassOf']/@*[local-name()='resource']" doid.owl | cut -d\" -f2`} />

Now the output only contains the URIs.

We can now create a script that receives multiple URIs given as standard input and the OWL file where to find all the parents as argument. The script named getparents.sh

<Execute command="nano getparents.sh" />

should contain the following lines:

```bash
OWLFILE=$1
xargs -I {} xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='{}']/*[local-name()='subClassOf']/@*[local-name()='resource']" $OWLFILE | \
cut -d\" -f2
```

and add the permission:

<Execute command="chmod u+x getparents.sh" />

To get the parents of malignant hyperthermia, we will only need to give the URI as input and the OWL file as argument:

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_8545' | ./getparents.sh doid.owl" />

The output will include the URIs of the two parents.

#### Labels of parents

But if we need the labels we can redirect the output to the `getlabels.sh` script:

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_8545' | ./getparents.sh doid.owl | ./getlabels.sh doid.owl" />

The output will now be the label of the parents of malignant hyperthermia.

Again, the same can be done with caffeine:

<Execute command="echo 'http://purl.obolibrary.org/obo/CHEBI_27732' | ./getparents.sh chebi_lite.owl | ./getlabels.sh chebi_lite.owl" />

And now the output contains the labels of the parents of caffeine.

#### Related classes

If we are interested in using all the related classes besides the ones that represent a generalization (`subClassOf` ), we have to change our XPath to:

<Execute command={`xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/CHEBI_27732']/*[local-name()='subClassOf']//*[local-name()='someValuesFrom']/@*[local-name()='resource']" chebi_lite.owl | cut -d\" -f2`} />

We should note that these related classes are in the attribute resource of `someValuesFrom` element inside a subClassOf element. The URIs of the related classes of caffeine are now displayed.

#### Labels of related classes

To get the labels of these related classes, we only need to add the getlabels.sh
script:

<Execute command={`xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/CHEBI_27732']/*[local-name()='subClassOf']//*[local-name()='someValuesFrom']/@*[local-name()='resource']" chebi_lite.owl | cut -d\" -f2 | ./getlabels.sh chebi_lite.owl`} />

The output is now terms that we could use to expand our text processing.

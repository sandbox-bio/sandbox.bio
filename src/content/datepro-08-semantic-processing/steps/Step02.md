<script>
import Execute from "$components/Execute.svelte";
</script>

##  URIs and Labels

In the previous examples, we searched the OWL file using labels and URIs. To standardize the process, we will create two scripts that will convert a label into a URI and vice-versa. The idea is to perform all the internal ontology processing using the URIs and in the end convert them to labels, so we can use them in text processing.

#### URI of a label

To get the URI of malignant hyperthermia, we can use the following query:

<Execute command="xmllint --xpath "//*[local-name()='label' and text()='malignant hyperthermia']/../@*[local-name()='about']" doid.owl" />

We added the `@*[local-name()='about']` to extract the URI specified as an attribute of that class.
The output will be the name of the attribute and its value:
```xml
rdf:about="http://purl.obolibrary.org/obo/DOID_8545"
```

To extract only the value, we can add the string function to the XPath query:

<Execute command="xmllint --xpath "string(//*[local-name()='label' and text()='malignant hyperthermia']/../@*[local-name() ='about'])" doid.owl" />

The output will now be only the attribute value.
Unfortunately, the string function returns only one attribute value, even if many are matched. Nonetheless, we use the string function because we assume that malignant hyperthermia is an unambiguous label, i.e. only one
class will match. To avoid this limitation we can add the cut command using the character delimiting the URI, i.e. ".

<Execute command="xmllint --xpath "//*[local-name()='label' and text()='malignant hyperthermia']/../@*[local-name()='about']" doid.owl | cut -d\" -f2" />

To get the URI of caffeine is just about the same command:

<Execute command="xmllint --xpath "//*[local-name()='label' and text()='caffeine']/../@*[local-name()='about']" chebi_lite.owl | cut -d\" -f2" />

We can now write a script that receives multiple labels given as standard input and the OWL file where to find the URIs as argument. Thus, we can create the script named `geturi.sh` 

<Execute command="nano geturi.sh" />

with the following lines:
```bash
OWLFILE=$1
xargs -I {} xmllint --xpath "//*[local-name()='label'
and text()='{}']/../@*[local-name()='about']" $OWLFILE | \
cut -d\" -f2
```

Again we cannot forget to save the file in our working directory, and add the right permissions using `chmod` as we did with our scripts in the previous turorials. 

<Execute command="chmod u+x geturi.sh" />

The `xargs` command is used to process each line of the standard input.
Now to execute the script we only need to provide the labels as standard input:

<Execute command="echo 'malignant hyperthermia' | ./geturi.sh doid.owl" />

<Execute command="echo 'caffeine' | ./geturi.sh chebi_lite.owl" />

The output should be the URIs of those classes.

We can also execute the script using multiple labels, one per line:

<Execute command="echo -e 'malignant hyperthermia\nmuscle tissue disease' | ./geturi.sh doid.owl" />

<Execute command="echo -e 'caffeine\npurine alkaloid\ntrimethylxanthine' | ./geturi.sh chebi_lite.owl" />

The output will be a URI for each label.

#### Label of a URI

To get the label of the disease entry with the identifier `8545`, we can also use the `xmllint` command:

<Execute command="xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/DOID_8545']/*[local-name()='label']/text()" doid.owl" />

We added the `@*[local-name()='label']` to select the element within the class that describes the label.
The output should be the label we were expecting: malignant hyperthermia.

We can do the same to get the label of the compound entry with the identifier `27732`:

<Execute command="xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='http://purl.obolibrary.org/obo/CHEBI_27732']/*[local-name()='label']/text()" chebi_lite.owl" />

Again, the output should be the label we were expecting:
caffeine. We can now write a script that receives multiple URIs given as standard input and the OWL file where to find the labels. We can create a script named getlabels.sh 

<Execute command="nano getlabels.sh" />

with the following lines:

```bash
OWLFILE=$1
xargs -I {} xmllint --xpath "//*[local-name()='Class'][@*[local-name()='about']='{}']/*[local-name()='label']/text()" $OWLFILE
```

The `xargs` command is used to process each line of the standard input.  Now to execute the script we only need to provide the URIs as standard input:

<Execute command="chmod u+x getlabels.sh" />

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_8545' | ./getlabels.sh doid.owl" />

<Execute command="echo 'http://purl.obolibrary.org/obo/CHEBI_27732' | ./getlabels.sh chebi_lite.owl" />

The output should be the labels of those classes:
- malignant hyperthermia
- caffeine

We can also execute the script with multiple URIs:

<Execute command="echo -e 'http://purl.obolibrary.org/obo/DOID_8545\nhttp://purl.obolibrary.org/obo/DOID_66' | ./getlabels.sh doid.owl" />

<Execute command="echo -e 'http://purl.obolibrary.org/obo/CHEBI_27732\nhttp://purl.obolibrary.org/obo/CHEBI_26385\nhttp://purl.obolibrary.org/obo/CHEBI_27134' | ./getlabels.sh chebi_lite.owl" />


The output will be a label for each URI.

To test both scripts, we can feed the output of one as the input of the other, for example:

<Execute command="echo -e 'malignant hyperthermia\nmuscle tissue disease' | ./geturi.sh doid.owl | ./getlabels.sh doid.owl" />

<Execute command="echo -e 'caffeine\npurine alkaloid\ntrimethylxanthine' | ./geturi.sh chebi_lite.owl | ./getlabels.sh chebi_lite.owl" />

The output will be the original input, i.e. the labels given as arguments to the `echo` command.

Now we can use the URIs as input:

<Execute command="echo -e 'http://purl.obolibrary.org/obo/DOID_8545\nhttp://purl.obolibrary.org/obo/DOID_66' | ./getlabels.sh doid.owl | ./geturi.sh doid.owl" />

<Execute command="echo -e 'http://purl.obolibrary.org/obo/CHEBI_27732\nhttp://purl.obolibrary.org/obo/CHEBI_26385\nhttp://purl.obolibrary.org/obo/CHEBI_27134' | ./getlabels.sh chebi_lite.owl | ./geturi.sh chebi_lite.owl" />

Again the output will be the original input, i.e. the URIs given as arguments to the `echo` command.
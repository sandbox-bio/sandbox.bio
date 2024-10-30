<script>
import Execute from "$components/Execute.svelte";
</script>

#### Complex Elements

Not always the XML elements are in the same line, as fortunately was the case
of the PubMed identifiers. In those cases, we may have to use the `xmllint`
command, a parser that is able to extract data through the specification of a
XPath query, instead of using a single line pattern as in `grep`.

#### XPath

XPath (XML Path Language) is a powerful tool to extract information from
XML and HTML documents by following their hierarchical structure. Check
W3C for more about XPath syntax.

##### Namespace problems

In the case of our protein XML files, we can see that their second line defines
a specific namespace using the `xmlns` attribute:

```xml
<uniprot
    xmlns="http://uniprot.org/uniprot"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://uniprot.org/uniprot
    http://www.uniprot.org/support/docs/uniprot.xsd">
```

This complicates our XPath queries, since we need to explicitly specify that
we are using the local name for every element in a XPath query. For example,
to get the data in each reference element:

<Execute command={`xmllint --nsclean --xpath "//*[local-name()='reference']" chebi_27732_P21817.xml`} />

We should note that `//` means any path in the XML file until reaching a
reference element. The square brackets in XPath queries normally represent
conditions that need to be verified. The `--nsclean` removes the redundant
namespace declaration in entry.

#### Only local names

If we are only interested in using local names there is a way to avoid the
usage of local-name() for every element in a XPath query. We can identify
the top-level element, in our case entry, and extract all the data that it
encloses using a XPath query. For example, we can create the auxiliary file
`chebi_27732_P21817_entry.xml` by adding the redirection operator:

<Execute command={`xmllint --nsclean --xpath "//*[local-name()='entry']" chebi_27732_P21817.xml > chebi_27732_P21817_entry.xml`} />

The new XML file now starts and ends with the entry element without
any namespace definition:

```xml
<entry dataset="Swiss-Prot" created="1991-05-01" ...
<accession>P21817</accession>
...
</sequence>
</entry>
```

Now we can apply any XPath query, for example `//reference`, on the
auxiliary file without the need to explicitly say that it represents a local name:

<Execute command={`xmllint --xpath '//reference' chebi_27732_P21817_entry.xml`} />

The output should contain only the data inside of each reference element.

#### Queries

The XPath syntax allow us to create many useful queries, such as:

- `//dbReference` - elements of type dbReference that are descendants
  of something:

<Execute command={`xmllint --xpath '//dbReference' chebi_27732_P21817_entry.xml`} />

- `/entry//dbReference` - equivalent to the previous query but specifying that the dbReference elements are descendants of the entry element:

<Execute command={`xmllint --xpath '/entry//dbReference' chebi_27732_P21817_entry.xml`} />

- `/entry/reference/citation/dbReference` - similar to the previous query but specifying the full path in the XML file, i.e. only `dbReference`
  elements descendants of citation, reference and entry elements:

<Execute command={`xmllint --xpath '/entry/reference/citation/dbReference' chebi_27732_P21817_entry.xml`} />

- `//dbReference/*` - any child elements of a dbReference element:

<Execute command={`xmllint --xpath '//dbReference/*' chebi_27732_P21817_entry.xml`} />

- `//dbReference/property[1]` - first property element of each dbReference
  element:

<Execute command={`xmllint --xpath '//dbReference/property[1]' chebi_27732_P21817_entry.xml`} />

- `//dbReference/property[2]` - second property element of each
  dbReference element:

<Execute command={`xmllint --xpath '//dbReference/property[2]' chebi_27732_P21817_entry.xml`} />

- `//dbReference/property[3]` - third property element of each dbReference
  element:

<Execute command={`xmllint --xpath '//dbReference/property[3]' chebi_27732_P21817_entry.xml`} />

- `//dbReference/property/@type` - all type attributes of the property
  elements:

<Execute command={`xmllint --xpath '//dbReference/property/@type' chebi_27732_P21817_entry.xml`} />

- `//dbReference/property[@type="protein sequence ID"]` - the
  previous property elements that have an attribute type equal to protein
  sequence ID:

<Execute command={`xmllint --xpath '//dbReference/property[@type="protein sequence ID"]' chebi_27732_P21817_entry.xml`} />

- `//dbReference/property[@type="protein sequence ID"]/@value` - the string assigned to each attribute value of the previous property
  elements:

<Execute command={`xmllint --xpath '//dbReference/property[@type="protein sequence ID"]/@value' chebi_27732_P21817_entry.xml`} />

- `/entry/sequence/text()` - the contents inside the sequence element:

<Execute command={`xmllint --xpath '/entry/sequence/text()' chebi_27732_P21817_entry.xml`} />

Thus, an alternative way to extract the PubMed identifiers using `xmllint`
instead of `grep`, would be something like this:

<Execute command={`xmllint --xpath '//dbReference[@type="PubMed"]/@id' chebi_27732_P21817_entry.xml`} />

However, the output contains all identifiers in the same line and with the
id label.

Previous versions of xmllint may print all the output in the same line.
In that case, we need to add an extra `tr ' ' '\n'` command to split the
output in multiple lines (one line per identifier).
Extracting XPath results
To extract the identifiers, we can use the cut command:

Extracting XPath results
To extract the identifiers, we can use the cut command:

<Execute command={`xmllint --xpath '//dbReference[@type="PubMed"]/@id' chebi_27732_P21817_entry.xml | cut -d'"' -f2`} />

The cut command extracts the value inside the double quotes.

<script>
import Alert from "$components/Alert.svelte";
</script>

This is the **fourth tutorial** in a series that will demonstrate how shell scripting can be used to perform the tasks that health and life science specialists may need to undertake to find and retrieve biomedical data and text. We will use the compound caffeine as an example and explore different public repositories to identify diseases related to it. The focus is not on the specific relationships we may discover, but on the process of obtaining them.

The objective of this tutorial is to learn how to use task repetition techniques to efficiently apply the same task to all proteins in the list we gathered in the previous tutorial. The task will be to construct the URL that will enable us to download the information from UniProt for each protein using the [RESTful web services provided by UniProt](https://www.uniprot.org/help/api).

<Alert>
This tutorial is part of a series of tutorials adapted as interactive versions of the hands-on steps described in the [Data and Text Processing for Health and Life Sciences](https://labs.rd.ciencias.ulisboa.pt/book/) book, which is licensed under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). 
</Alert>
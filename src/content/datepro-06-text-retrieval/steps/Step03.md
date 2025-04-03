<script>
import Execute from "$components/Execute.svelte";
</script>

Disease Recognition

Instead of reading all that text to find any disease related with caffeine, we
can try to find sentences about a given disease by using grep:

<Execute command="grep 'malignant hyperthermia' chebi_27732.txt" />

To save the filtered text in a file named `chebi_27732_hyperthermia.txt`, we
only need to add the redirection operator:

<Execute command="grep 'malignant hyperthermia' chebi_27732.txt > chebi_27732_hyperthermia.txt" />

This is a very simple way of recognizing a disease in text. The next tutorials
will describe how to perform more complex text processing tasks.

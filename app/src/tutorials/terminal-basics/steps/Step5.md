<script>
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

We can store the result of our analyses in files using the `>` operator.

For example, to store all chicken burrito orders into `burritos.tsv`:

<Execute command='grep "Chicken Burrito" orders.tsv > burritos.tsv' />

Now notice that <Execute command="ls" inline /> shows the newly-created file!

---

<Alert color="warning">

Be careful when using `>`. If you have a FASTA file and you want to extract lines with the character `>`, make sure to put it in quotes! Otherwise, grep interprets `>` as an operator and will truncate your FASTA file.

This is a rite of passage in bioinformatics but can be avoided if you prepend `cat` to your pipelines, e.g. `cat ref.fa | grep ">"` (this is why there is no such thing as a _useless use of cat_).

<!-- Let's illustrate this pitfall in our sandbox. If we grep for `>` with quotes, it works as expected:
<Execute command='grep ">" copy.fa' />
Whereas not quoting results in a truncated file!
<Execute command='grep > copy.fa' />
You can verify that using `ls`:
<Execute command='ls' /> -->
</Alert>

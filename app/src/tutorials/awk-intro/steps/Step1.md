<script>
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

The first thing to know about `awk` is that it operates on **rows and columns**, also known as **records and fields**, respectively.

To extract the 3rd column of `orders.tsv`, we can use the `print` command and use `$n` to represent the `n`th column:

<Execute command={`awk '{ print $3 }' orders.tsv | head`} />

Note that `print $0` (or just `print`) refers to the entire line:

<Execute command={`awk '{ print $0 }' orders.tsv | head`} />

`Awk` also provides built-in variables: for example, `NF` (number of fields) tells you how many columns there are on the current line. In turn, `$NF` means fetch the last column of the file:

<Execute command={`awk '{ print $NF }' orders.tsv | head`} />

Let's now extract columns 2 and 3 to obtain the items ordered and the quantity for each:

<Execute command={`awk '{ print $2, $3 }' orders.tsv | head`} />

Note, however, that when you compare the output of the command above to <Execute command="head orders.tsv" inline />, it seems the `item_name` column is truncated! This is because by default, `awk` treats tabs and spaces both as column delimiters, so we have to explicitely tell it we have tab-separated data:

<Execute command={`awk -F "\\t" '{ print $2, $3 }' orders.tsv | head`} />

<Alert>
	Remember to explicitely tell `awk` what your delimiters are!
</Alert>

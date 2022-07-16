<script>
import Execute from "components/Execute.svelte";
</script>

So far, we only used variables inside the awk script. But what if we wanted to pass a Bash variable into an awk script? This can be useful if you want to write a general script that you can call with different parameters, or because your script depends on system parameters.

For example, let's define a variable that represents the tax rate:

<Execute command={`TAX_RATE=0.15`} />

We can now use `-v` to pass this variable into `awk`:

<Execute command={`awk -F "\\t" -v tax=$TAX_RATE ' \\ BEGIN { print tax }'`} />

To add the tax onto the price column (5th column), we first need to remove the dollar sign from the price so we can do arithmetic on it. Let's use `substr` to fetch everything after the first character:

<Execute command={`awk -F "\\t" -v tax=$TAX_RATE ' \\ { if(NR > 1) print $3, substr($5, 2) }' orders.tsv | head`} />

Now combine that with the variable to calculate the price including tax:

<Execute command={`awk -F "\\t" -v tax=$TAX_RATE ' \\ { if(NR > 1) print $3, substr($5, 2) * (1 + tax) }' orders.tsv | head`} />

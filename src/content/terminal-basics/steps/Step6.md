<script>
import Execute from "$components/Execute.svelte";
</script>

Finally, let's explore environment variables. For example, to define a variable `abc` with contents `123`, and `def` with `hello`:

<Execute command="abc=123" />

<Execute command="def=hello" />

Make sure there are no spaces surrounding the equal sign, otherwise the terminal treats the variable name as a command (try it!).

To display the content of the variable `abc`, use the `$` delimiter:

<Execute command="echo $abc" />

To delete a variable, use `unset`:

<Execute command="unset abc" />

Use `env` to view all available variables:

<Execute command="env" />

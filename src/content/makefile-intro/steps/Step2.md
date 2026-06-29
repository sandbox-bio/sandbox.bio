<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Makefile prefixes -->

You probably noticed that our example prints out the command before executing it, which
is not always ideal. Luckily, `make` gives us a couple of ways to manage how each
command is executed:

1. **`@`**: Tells `make` not to print the command before executing.
2. **`-`**: Tells `make` to ignore any errors and continue executing the recipe.

Change `echo` to `@echo` and rerun to see the difference.

<Execute command="make hello-world" />

> **Optional: Test `-`**  
> Add a second command to the `hello-world` recipe above the `echo` command and rerun:
>
> ```Makefile
> 	rm this-file-does-not-exist.txt
> ```
>
> Add the `-` prefix to `rm` and notice the difference when you re-run. After, be sure
> to remove the `rm` line.

<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Special variables -->

So far, we've told `make` which rule to run, but remember -- we loathe repeating
ourselves and typing more than necessary. It turns out `make` doesn't require us to
provide a name. Try and see what happens.

<Execute command="make" />

Without a target specified, `make` will run the first rule it sees... which is not
ideal. We can set a special variable to tell `make` which rule to run by default, but
first, we need to understand make variables.

Add the following to your file, then test it.

```Makefile
NAME := value

$(NAME):
	echo "$(NAME)"
```

<Execute command="make value" />

Essentially:

1. Declare variables at the top of the file using the pattern shown.
2. Use the variable anywhere in rules using `$(NAME)`.
3. The `:=` assignment operator evaluates the assignment immediately. There are other
   variations, but this one is all we'll need for now.

Special variables work just like regular variables, and `.DEFAULT_GOAL` sets the default
target to use if not specified. Add the following line to the top of the file and rerun.

```Makefile
.DEFAULT_GOAL := hello-world
```

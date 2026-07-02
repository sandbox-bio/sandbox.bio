<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Integrating the features -->

Go ahead and run `clean` then `story` to remind us where we stand:

<Execute command="make clean; make story" />

As expected, we get an error because the story files do not exist.

We can `make` the files with their individual targets, but `make` was _designed_ to
solve this exact problem and auto-run the rules if the targets don't exist!

We could explicitly declare the prerequisites for each phony rule, e.g.:

```Makefile
story-beginning: story-beginning.txt
```

Not ideal, but we can't use implicit rules with phony targets. Luckily, there's another
pattern we _can_ use.

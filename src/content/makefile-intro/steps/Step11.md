<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Pattern substitution -->

Unfortunately, the explicit vs implicit pattern we ran into has another caveat: `make`
will not evaluate implicit rules for phony targets. No implicit rules _yet_, but we
_can_ still streamline our rules using builtin functions.

Let's build some lists that define the target names we want to use. Add the following at
the top of `Makefile.story`:

```Makefile
STORY_FILES := $(wildcard story-*.txt)
STORY_PARTS := $(patsubst %.txt,%,$(STORY_FILES))
```

> _Make sure you don't include spaces between the commas -- those will be added to the
> pattern if you do!_

Basically, the above statements do the following:

1. Use `*` and `wildcard` to find all of the files that match the pattern.
2. Remove the `.txt` extension from each of the files we found.

We can then use these to replace each of the individual `story-` targets with the
variable:

```Makefile
.PHONY: story $(STORY_PARTS)
$(STORY_PARTS):
	@cat "${@}.txt"
```

Try running the updated makefile, and we'll dive into implicit rules next.

<Execute command="make" />

> _If your change didn't print the story as expected, did you remember to keep the rules
> that define the story order prerequisites?_
>
> ```Makefile
> story-middle: story-beginning
> story-end: story-middle
> story: story-end
> 	@echo "The End."
> ```

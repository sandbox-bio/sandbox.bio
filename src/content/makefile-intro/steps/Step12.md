<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Variables -->

In order to use implicit rules, we need file targets. We could make individual file
targets like this -- using the special variable `$@` in place of rewriting the target
name:

```Makefile
story-beginning.txt:
	@echo "Once upon a time, there was a chicken." > $@
```

But we'd have to write a separate rule for each, and that just sounds tedious. Instead,
we can define the story in variables:

```Makefile
BEGINNING := "Once upon a time, there was a chicken."
MIDDLE := "It crossed the road."
END := "And lived happily ever after."
```

And then use a shorthand for `patsubst` called substitution reference to define the
target names:

```Makefile
STORY_PARTS := beginning middle end
STORY_FILES := $(STORY_PARTS:%=story-%.txt)
STORY_RULES := $(STORY_PARTS:%=story-%)
```

Let's add a sanity check rule and run it before we hook this into our makefile:

```Makefile
check:
	@echo "parts: ${STORY_PARTS}"
	@echo "files: ${STORY_FILES}"
	@echo "rules: ${STORY_RULES}"
```

<Execute command="make check" />

You should see the following printed out:

```shell
parts: beginning middle end
files: story-beginning.txt story-middle.txt story-end.txt
rules: story-beginning story-middle story-end
```

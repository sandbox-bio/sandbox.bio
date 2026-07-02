<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Static pattern rules -->

Luckily, `make` supports a feature called _static pattern rules_ that allow us to define
a single rule to handle the dependencies matching. These rules use the following
pattern:

```Makefile
targets ...: target-pattern: prerequisites ...
	recipe
```

The `target-pattern` works like the pattern matching with `%` we've used before. Find
your `$(STORY_PARTS)` rule and update it as shown. Now build the story and watch it
work:

```Makefile
$(STORY_RULES): story-%: story-%.txt
	@cat "story-${*}.txt"
```

<Execute command="make clean; make story" />

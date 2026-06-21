<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Prerequisites Setup -->

We'll cover prerequisites over the next several steps.

We've already had enough greetings for one tutorial, so let's tell a story instead. And
what does every good story need? A beginning, middle, and end.

First, write out each part of the story into separate files.

<Execute command={`echo "Once upon a time, there was a chicken." > story-beginning.txt`}
/>

<Execute command={`echo "It crossed the road." > story-middle.txt`} />

<Execute command={`echo "And lived a happy life on the other side." > story-end.txt`} />

Switch to `Makefile.story` and add the following:

```Makefile
.PHONY: story story-beginning story-middle story-end
story-beginning:
	@cat "story-beginning.txt"

story-middle:
	@cat "story-middle.txt"

story-end:
	@cat "story-end.txt"

story: story-beginning story-middle story-end
	@echo "The End."
```

Notice:

1. **We're still using phony targets.** That's intentional, and we'll clean it up later.
2. **We declared `story` with three prerequisites**.

Try rerunning make with the new target:

<Execute command="make story" />

It failed! Let's go to the next step and fix it.

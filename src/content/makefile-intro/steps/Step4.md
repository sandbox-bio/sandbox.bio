<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Prerequisites -->

Let's try to define our own prequisites, but we've already had enough greetings for one
tutorial, so let's tell a story instead. And what does every good story need? A
beginning, middle, and end.

First, write out each part of the story.

<Execute command={`echo "Once upon a time, there was a chicken." > story-beginning.txt`}
/>

<Execute command={`echo "It crossed the road." > story-middle.txt`} />

<Execute command={`echo "And lived a happy life on the other side." > -storyend.txt`} />

Next, let's start a new Makefile and add a story:

<Execute command={`vim Makefile.story`} />

```Makefile
.DEFAULT_GOAL := story
.PHONY: story story-%
story-beginning:
	@cat "story-beginning.txt"

story-middle:
	@cat "story-middle.txt"

story-end:
	@cat "story-end.txt"

story: story-beginning story-middle story-end
	@echo "The End."
```

Two things to note before we run:

1. **We used the `.DEFAULT_GOAL` special variable.** The named target will be run
   automatically when calling `make` without a specific target.
2. **We used a wildcard to declare `story-%` as phony.** The `%` wildcard is used for
   string pattern matching in targets and prerequisites while the `*` wildcard can be
   used for file globbing.

<!-- prettier-ignore -->
<Quiz id="q4.1" choices={[
{ valid: false, value: "Only \"The End\""},
{ valid: false, value: "\"The End\" first, then the story lines."},
{ valid: false, value: "The story lines in order followed by \"The end\""},
{ valid: true, value: "The story lines randomly followed by \"The end\""},
]}>
<span slot="prompt">
What will be printed to stdout when we run `make`?
</span>
</Quiz>

We added a suffix to the filename, so we need to tell `make` which file to execute via
the `-f` flag.

<Execute command="make -f Makefile.story" />

Now, chances are good that the story was printed to stdout in order -- and that it will
continue to do so. However, there is _no_ guarantee that the prerequisites will be
called in the order they are defined.

We can easily guarantee the order by updating the `story-%` targets. Make the update and
re-execute `make` to see if anything changes.

```Makefile
story-middle: story-beginning
story-end: story-end
```

<Execute command="make -f Makefile.story" />

While we're at it, let's update the `story` target to use the pattern, too. No point in
being overly verbose if we don't need to, right? Update and run it again.

```Makefile
story: story-%
```

Not what you expected? It turns out that `make` differentiates between two types of
targets: concrete and dynamic. It boils down to wildcards in the target -- if we use a
wildcard in the target, it's a dynamic target, and we can use pattern matching in the
prerequisites. So far, we've only used concrete targets and special targets, but we'll
cover dynamic targets next.

## Recap

1. Use `.DEFAULT_GOAL := <TARGET>` to tell make which target to build if no targets are
   specified. Otherwise, the first target defined will be used
2. Use `%` for string matching and `*` for file name matching.
3. Prerequisites are not guaranteed to run in a specific order.

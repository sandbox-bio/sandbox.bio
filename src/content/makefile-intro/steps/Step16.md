<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Cleanup -->

In the last step we ran into the "TARGET is up to date" message again. While we could
use `-B` to force `make` to rebuild targets, we would likely end up building more rules
that we actually need or want -- not an issue for this example, but problematic if
you're compiling large libraries.

Instead, let's create a rule to clean up after ourselves.

<!-- prettier-ignore -->
<Quiz id="q16.1" choices={[
{ valid: false, value: "Mark it as the default group."},
{ valid: true, value: "Mark it as phony."},
{ valid: true, value: "Ignore errors with `-`."},
{ valid: false, value: "Suppress the command echo with `@`."},
{ valid: false, value: "Add prerequisites for the `story-%.text` targets"},
]}>
<span slot="prompt">
What will need for a `clean` rule?
</span>
</Quiz>

Because `clean` isn't a file, we'll need to mark it as phony, and because `rm` errors if
the file doesn't exist, we should ignore errors so the rest of the recipe can complete.

We'll also create a `clean-story` target: `clean` is a common target, and we are working
across two makefiles already. To ensure we don't clobber `clean` targets defined in our
other makefile, we'll list our specific target as a dependency of `clean`.

Add the following rule to `Makefile.story`, then build and check the file contents.

```Makefile
.PHONY: clean clean-story
clean: clean-story
clean-story:
	-@rm -f story-*.txt
```

<Execute command="make clean; make story-beginning.txt; cat story-beginning.txt" />

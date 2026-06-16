<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Dynamic Targets -->

Let's go back to the story makefile we were working on before.

<Execute command="cat Makefile.story" />

Do you hate retyping the same thing over and over, too? Thankfully, `make` supports a
few special variables -- one in particular that we can use to make this a little
cleaner: `$@`. This variable can be used within a recipe to refer to the name of the
target.

We can replace all of the `story-%` targets with a single, dynamic target. Don't forget
-- we _also_ need to keep the prerequisites to ensure the story order.

<Execute command="vim Makefile.story" />

```Makefile
story-%:
	@cat "${@}.txt"

story-middle: story-beginning
story-end: story-middle
```

Before we run it, let's update our original `Makefile` so that we don't have to declare
`-f Makefile.story` every time we call `make`.

<Execute command="vim Makefile" />

```Makefile
include Makefile.story
```

Now, when we call make, we should see our short story printed out.

<Execute command="make" />

## Recap

- Use wildcards like `%` and special variables like `$@` to reduce repitition within
  targets and recipes.
- Include one makefile within another via the `include` directive at the top of your
  file.

<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Include -->

We created a new `Makefile.story`, but `make story` failed. Why?

Although `make` will work with any properly formatted file, it expects `Makefile` by
default. We can tell `make` to use a different file via the `-f` flag:

<Execute command="make -f Makefile.story" />

<!-- prettier-ignore -->
<Quiz id="q7.1" choices={[
{ valid: false, value: "It will print out the full story."},
{ valid: true, value: "It will print out one part of the story"},
{ valid: false, value: "It will error."},
{ valid: false, value: "Nothing to do for the targets"},
]}>
<span slot="prompt">
What do you think will happen when we run the above command?
</span>
</Quiz>

Running the above command should execute the first target in the file -- the beginning.
But having to type `-f Makefile.story` every time is tiresome. Do we really need to type
this every time we run `make`? Nope! There's a better solution: Tell `make` to include
our new file.

Go back to `Makefile` and make the following edits:

1. Add the following at the top of the file (before `.DEFAULT_GOAL`)

   ```Makefile
   include Makefile.story
   ```

2. Update `.DEFAULT_GOAL` to use `story` as the default target, too.

`make` will complain if multiple included files define a `DEFAULT_GOAL`, so it's best to
include files before we define the default, which ensures the value in the primary file
overrides all the others.

Now, we can run `make` as before to build the target:

<Execute command="make" />

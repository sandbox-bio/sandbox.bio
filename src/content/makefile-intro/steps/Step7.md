<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Include -->

We created a new `Makefile.story`, but `make story` failed. Why?

Although `make` will work with any properly formatted file, it expects `Makefile` by
default. We can tell `make` to use a different file via the `-f` flag:

<Execute command="make -f Makefile.story" />

But once again, do we really need to type this every time we run `make`? Nope! There's a
better solution: Tell `make` to include our new file.

Go back to `Makefile` and make the following edits:

1. Add the following at the top of the file (before `.DEFAULT_GOAL`)

   ```Makefile
   include Makefile.story
   ```

2. Update `.DEFAULT_GOAL` to use `story` as the default target.

`make` will complain if multiple included files define a `DEFAULT_GOAL`, so it's best to
include files before we define the default, which ensures the value in the primary file
overrides all the others.

Now, we can run `make` as before to build the target:

<Execute command="make" />

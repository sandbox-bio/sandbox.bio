<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: PHONY -->

Thankfully, there are two ways to force `make` to execute a rule. The first is to simply
force make to run the target with `-B`:

<Execute command="make -B hello-world" />

But we don't like to type extra flags every time... and we don't want to force every
rule to build every time. Instead, we can declare the target as "fake" so that `make`
will skip the file check using the special target `.PHONY`.

Add the following to the line above the `hello-world` rule and re-run the command

```Makefile
.PHONY: hello-world
```

<Execute command="make hello-world" />

`.PHONY` is just one of the special targets supported by `make`, and we added our target
as a prerequisite. We'll cover more special targets and prerequisites later, but a
couple of things to note now:

1. It doesn't matter where you put `.PHONY`, but the rule of thumb is to put it at the
   top of the file or just before the rule(s) it affects.
2. You can declare multiple `.PHONY` targets. In fact, you can decalre multiple of any
   target: `make` will merge the prerequisites of each, but it will only execute the
   _last_ recipe.

<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Makefile structure. Cover targets and recipes. -->

A `Makefile` defines a series of rules with the general format:

```Makefile
target: prerequisites
	recipe
```

1. **Target**: The name of the rule.
2. **Prerequisites**: (Optional) The names of other rules that this rule depends on.
3. **Recipe**: (Optional) One or more shell commands that tell `make` how to build the
   target.
   1. MUST have a single tab indent (no spaces)

Take a look at the `Makefile`, which contains the quintessential example: The universal
"Hello, world". Tell make to build the rule and see what happens:

<Execute command="make hello-world" />

You should see "Hello, world" printed to the console. Congrats! You've made your first
Makefile.

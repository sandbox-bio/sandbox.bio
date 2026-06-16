<!-- <script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Makefile structure. Cover targets and recipes. -->

> **VIM Shortcuts**
>
> - `i`: INSERT mode
> - `:wq`: Save and close the editor

Let's start with a simple example: The universal "Hello, world".

> NOTE: Sometimes, the `Makefile` is created with a full HTML page, so we'll clear it
> out first.

<Execute command={`echo '' > Makefile; vim Makefile`} />

Write the following, then save and close the editor.

```Makefile
hello-world:
	echo -e "Hello, world\\n"
```

Let's breakdown what you just wrote:

1. **The target**; i.e., text at the start of a line followed by a colon; e.g.,
   `hello-world`. This assigns a _name_ to a recipe, so that you can tell `make` which
   recipe to, well, _make_.
2. **The recipe**; i.e., a line starting with a tab (`\t`) followed by a shell command;
   e.g., `echo ...`. This tells `make` what to do for a given target.

Some important details:

- A target MUST start at the beginning of the line and MUST be followed by a colon.
- A recipe MUST follow a target and MUST start with a single tab character (no spaces as
  tab).
- A target does not need a recipe, and a recipe may be multiple lines long.

Now execute the target.

<Execute command="make hello-world" />

You should see "Hello, world" printed to the console. Congrats! You've made your first
Makefile.

You also probably saw `echo "Hello, world"`: `make` will print each command before
running it. Usually, this is great, but in some caeses, like with `echo`, it's tedious.

Fortunately, there's a simple fix: Put an `@` just before the command.

```Makefile
	@echo ...
```

#### Recap

- Define **targets** as the name of the thing to run.
- Define **recipes** as a series of commands to run when executing a target.
- Use an `@` at the start of recipe lines to avoid printing the command to stdout. -->

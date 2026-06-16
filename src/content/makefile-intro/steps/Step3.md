<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: PHONY vs file targets -->

> **VIM Shortcuts**
>
> - `i`: INSERT mode
> - `:wq`: Save and close the editor

There's a small issue with the `hello-world` target we added in the previous step. Let's
break it!

<Execute command={`echo "Goodbye, world" > hello-world`} />

We created a file with the same name as our target.

<!-- prettier-ignore -->
<Quiz id="q3.1" choices={[
{ valid: true, value: "It will error or warn."},
{ valid: false, value: "It will print \"Hello, world\""},
{ valid: false, value: "It will print \"Goodbye, world\""},
{ valid: false, value: "Nothing"},
]}>
<span slot="prompt">
What do you think will happen when we run the make target?
</span>
</Quiz>

Let's see what happens if we try run the target now.

<Execute command="make hello-world" />

Originally, `make` was developed to automate code compiling in C. The goal was to
automatically detect files that had changes and rebuild the relevant binaries. By
default, `make` assumes that every _target_ is actually a _file_.

In other words, the name "make" is quite literal: The goal is to "make" the file
"hello-world". When we created the file, `make` reads the target, sees the file, and
realizes it has nothing to do. Great for automation, but it would be nice to guarantee a
target will always run.

Thankfully, there are two ways around this. The first is to simply force make to run the
target with `-B`

<Execute command="make -B hello-world" />

But we don't want to include `-B` every time -- that's too many letters (pun intended).
Instead, we can declare the target as "fake" so that `make` will skip the file check.

Let's edit the Makefile.

<Execute command="vim Makefile" />

On the line before the target, add the following:

```Makefile
.PHONY: hello-world
```

Re-run the make command, and you should see the "Hello, world" message show up again.

<Execute command="make hello-world" />

There are two important things going on here:

1. First, we declared a new target `.PHONY`.

   It turns out that `make` defines several special targets that allow us to customize
   and manage behavior. `.PHONY` simply tells `make` that these targets don't actually
   exist and to run them unconditionally. There are several others, but we won't cover
   them here. Google "makefile special targets" if you want to learn more.

2. Second, we declared `hello-world` as a prerequisite.

   Prerequisites are the key feature that make `make` so powerful. They tell `make` to
   execute them before the target.

#### Recap

1. `make` assumes that every target refers to an actual file.
2. `make` checks to see if a target is up to date before executing.
3. We can use the `.PHONY` special target to tell `make` that certain targets are not
   actually files and should always be run.

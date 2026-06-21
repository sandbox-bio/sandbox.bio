<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: DIY functions -->

In earlier versions of make, converting to upper or lower case is fairly easily done,
but it requires a few steps. Go ahead and add the solution below to the top of your
Makefile, and we'll walk through each part.

```Makefile
to_upper = $(shell echo "$(1)" | tr '[:lower:]' '[:upper:]')
```

###### `to_upper =`

- We're essentially defining a new function we can use. Note the use of `=` instead of
  `:=`. While the latter evaluates immediately, the former evaluates every time it is
  used.

###### `$(shell ... | ...)`

- `shell` is a special function that allows us to execute shell commands within a
  variable (or more specifically, outside of a recipe). Everything after is treated as
  it would be on the command line.
- The pipe `|` character behaves exactly as it does on the command line: Separates two
  commands and passes the output of one as input to the other.

###### `echo "$(1)"`

- Print the first input parameter `$(1)`

###### `tr '[:lower:]' '[:upper:]'`

- `tr` (_translate_) is a bash command that processes input character by character
- Here, we're converting every lower case character into upper case.

Now, we can update our implicit rule to use `to_upper` with the special function `call`,
then build the recipe again and print out the file:

```Makefile
story-%.txt:
	@echo "$($(call to_upper,$*))" > $@
```

<Execute command="make story-beginning.txt; cat story-beginning.txt" />

Did it print "\* is up to date"? Let's make it easier to reset our targets next.

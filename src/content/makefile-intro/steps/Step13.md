<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: File targets -->

In the last step we added variables for the story sections and targets. Let's add rules
to process them.

We want rules for each `*.txt` file that will write the story section to the
appropriately named file. Should be pretty easy, right?

- We know we want to use implicit rules, so `story-%.txt:` for the targets.
- We can use `echo ... > $@` to write to the file based on the target name.
- We stored the story lines in variables, but can we use the target to specify which
  variable to use?

Yes! We can! The special vairable `$*` is the value of the pattern matched by `%`, and
we can wrap that variable in `$( )`, and `make` will evaluate the inner expression
first!

```Makefile
story-%.txt:
	@echo "$($*)" > $@
```

Save and execute one of the story targets, then `cat` the file out.

<Execute command="make story-beginning.txt; cat story-beginning.tx" />

However, there's a catch. Did you see it? It creates the file, but it's empty!

We used upper case for the variable name and lower case for the filename, so `$($*)`
will evaluate to `end` instead of `END`. If you're using `make` v4.4, there's a builtin
function for you: `$(upper VALUE)`, but we're going to do it the hard way.

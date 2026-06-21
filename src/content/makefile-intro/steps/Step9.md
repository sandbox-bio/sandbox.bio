<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Wildcards -->

When we created `Makefile.story`, you may have noticed the use of `story-%` within the
phony target. `make` supports a couple of wildcard options:

- `%`: Expands for string matching (typically, use this)
- `*`: Expands for file globbing

We can use these to write **implicit** rules that trigger based on a pattern instead of
the **explicit** rules we've been using so far. In the last step, we ensured the rule
order with prerequisites, so let's try out wildcards in the prerequisites for `story`:

```Makefile
story: story-%
```

<!-- prettier-ignore -->
<Quiz id="q3.1" choices={[
{ valid: false, value: "It will error or warn."},
{ valid: false, value: "It will print the whole story."},
{ valid: true, value: "It will print \"The End\""},
{ valid: false, value: "Nothing"},
]}>
<span slot="prompt">
What do you think will happen when we run the make target?
</span>
</Quiz>

<Execute command="make" />

Did it do what you expected? We'll discuss what happend and the fix in the next step.

<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Wildcards -->

As I'm sure you noticed, there is fair amount of redundancy in `Makefile.story` with all
of the `story-` names. With larger projects, this redundancy would lower readability and
maintainability. Fortunately, `make` supports a couple of wildcard options:

- `%`: Expands for string matching (typically, use this)
- `*`: Expands for file globbing

We can use these to write **implicit** rules that trigger based on a pattern instead of
the **explicit** rules we've been using so far. In the last step, we ensured the rule
order with prerequisites, so let's try out wildcards in the prerequisites for `story`:

```Makefile
story: story-%
```

<!-- prettier-ignore -->
<Quiz id="q9.1" choices={[
{ valid: true, value: "It will error or warn."},
{ valid: false, value: "It will print the whole story."},
{ valid: false, value: "It will print \"The End\""},
{ valid: false, value: "Nothing"},
]}>
<span slot="prompt">
What do you think will happen when we run the make target?
</span>
</Quiz>

<Execute command="make" />

Did it do what you expected? We'll discuss what happened and how to fix it in the next
step.

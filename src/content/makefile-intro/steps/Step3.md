<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

<!-- TOPIC: Make expects file targets -->
<!-- Rules are made to be broken -->

There's a small issue with our `hello-world` rule. Let's break it!

Create a file with the same name as our target.

<Execute command={`echo "Goodbye, world" > hello-world`} />

<!-- prettier-ignore -->
<Quiz id="q3.1" choices={[
{ valid: false, value: "It will error or warn."},
{ valid: false, value: "It will print \"Hello, world\""},
{ valid: false, value: "It will print \"Goodbye, world\""},
{ valid: true, value: "Nothing"},
]}>
<span slot="prompt">
What do you think will happen when we run the make target?
</span>
</Quiz>

Now try to run the target:

<Execute command="make hello-world" />

You should have gotten an information message saying that there was nothing to do for
the target. It turns out that `make` targets are expected to be files. In other words,
the name "make" is quite literal: The goal is to "make" the file "hello-world". When
`make` tries to build the rule, it sees the file and realizes it has nothing to do.
We'll fix this in the next step.

<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Dynamic vs. Concrete Targets -->

In the last step, we used wildcards to simplify the `story` prerequisites to `story-%`,
but `make` complained that it didn't have a rule for the `story-%` target. Remember what
we learned about **implicit** vs **explicit** rules? Well, it turns out `make`
differentiates between these types based on the _target_. Basically, if we use a
wildcard in the target, it's an implicit rule, and if not, it's an explicit rule.

**tl;dr:** `make` treated the `story-%` prerequisite as an explicit string, because the
`story` target is an explicit (static) target.

In this case, the fix is pretty simple: We already set prerequisites to ensure each rule
was executed in order, so `story` can depend only on `story-end`. Go ahead and change
`story-%` to `story-end` and check that the rule works again:

<Execute command="make" />

In the next step, we'll learn how do define a dynamic target to simplify our Makefile.

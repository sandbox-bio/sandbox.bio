<script>
import Execute from "$components/Execute.svelte";
</script>

<!-- TOPIC: Prereq Order -->

There's a subtle issue with our prerequisites: `make` does not guarantee that each rule
is run in the order defined.

Did you notice? Probably not. Chances are it'll run in the order we'd expect every time.

Still, if order _really_ matters, we can guarantee order by defining additional
prerequisites: end must depend on middle, and middle on beginning.

Update the targets as shown below and rerun to see if it works:

```Makefile
story-middle: story-beginning
story-end: story-end
```

<Execute command="make" />

Next, we'll cover patterns and implicit targets.

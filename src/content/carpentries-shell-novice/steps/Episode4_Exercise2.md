<script>
import Quiz from "$components/Quiz.svelte";
</script>

In our current directory, we want to find the 3 files which have the least number of
lines. Which command listed below would work?

<Quiz id="q4.2" choices={[
{ valid: false, value: "wc -l * > sort -n > head -n 3"},
{ valid: false, value: "wc -l * | sort -n | head -n 1-3"},
{ valid: false, value: "wc -l * | head -n 3 | sort -n"},
{ valid: true, value: "wc -l * | sort -n | head -n 3"},
]}>
<span slot="prompt">
</span>
</Quiz>

<script>
import Quiz from "$components/Quiz.svelte";
</script>

Suppose that you created a plain-text file in your current directory to contain a list of the statistical tests you will need to do to analyze your data, and named it `statstics.txt`. After creating and saving this file you realize you misspelled the filename!

<Quiz id="q3.2" choices={[
{ valid: false, value: "cp statstics.txt statistics.txt"},
{ valid: true, value: "mv statstics.txt statistics.txt"},
{ valid: false, value: "mv statstics.txt ."},
{ valid: false, value: "cp statstics.txt ."},
]}>
<span slot="prompt">
You want to correct the mistake, which of the following commands could you use to do so?
</span>
</Quiz>

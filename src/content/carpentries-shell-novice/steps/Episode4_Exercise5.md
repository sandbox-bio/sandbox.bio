<script>
import Quiz from "$components/Quiz.svelte";
</script>

Suppose you want to delete your processed data files, and only keep your raw files and processing script to save storage. The raw files end in .dat and the processed files end in .txt. Which of the following would remove all the processed data files, and only the processed data files?

<Quiz id="q4.5" choices={[
{ valid: false, value: "rm ?.txt"},
{ valid: true, value: "rm *.txt"},
{ valid: false, value: "rm * .txt"},
{ valid: false, value: "rm *.*"},
]}>
<span slot="prompt">
</span>
</Quiz>

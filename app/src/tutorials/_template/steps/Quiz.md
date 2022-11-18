<script>
import Quiz from "components/Quiz.svelte";
</script>

<!-- The "id" is used to maintain state on page refresh. It just has to be unique within a step in a tutorial -->
<Quiz id="q1" choices={[
	{ valid: true, value: "Shell"},
	{ valid: false, value: "Terminal"},
]}>
	<span slot="prompt">
		Is Bash a shell or a terminal?
	</span>
</Quiz>

<Quiz id="q2" choices={[
	{ valid: true, value: "ls -s -h Data"},
	{ valid: true, value: "ls -sh Data"},
	{ valid: false, value: "ls -size -h Data"},
	{ valid: true, value: "ls --size -h Data"},
	{ valid: false, value: "ls --sizeh Data"},
	{ valid: false, value: "ls --size-h Data"},
	{ valid: true, value: "ls -h -s Data"},
	{ valid: true, value: "ls -hs Data"},
	{ valid: false, value: "ls -hsize Data"},
]}>
	<span slot="prompt">
		Among the following commands, which ones are correct?
	</span>
</Quiz>

<Quiz id="q3" choices={[
	{ valid: true, value: "yes"},
	{ valid: false, value: "no"},
]}>
	<span slot="prompt">
		Now, type the following command in your terminal and then press <kbd>Enter</kbd> key: `date`

		Does the terminal display the current date?
	</span>
</Quiz>

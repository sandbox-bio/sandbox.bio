<script>
import Quiz from "components/Quiz.svelte";
</script>

The Shell is a software that handles Unix instructions. It is now commonly used by the bioinformatics community as it allows performing highly complex tasks by accessing either locally or distantly shared resources (e.g. files, databases, remote computing clusters...). These tasks can be coded in such a way that the multiple steps of a bioinformatics analysis can be saved into scripts. Those may then be shared with collaborators, traced and automated.

To access the Shell's powerful features, you need to learn a new way of thinking and a new language.

<Quiz id="q1" choices={[
	{ valid: true, value: "Yes"},
	{ valid: false, value: "No"},
]}>
	<span slot="prompt">
		Are you ready to learn Shell? (select "Yes" to continue)
	</span>
</Quiz>

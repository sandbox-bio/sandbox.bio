<script>
import Alert from "$components/Alert.svelte";
import Exercise from "$components/Exercise.svelte";
import Execute from "$components/Execute.svelte";

const criteria = [
{
	name: "File <code>chikungunya.jf</code> exists",
	checks: [{
		type: "file",
		path: "chikungunya.jf",
		action: "exists"
	}]
},
{
	name: "File <code>chikungunya.jf</code> contains the output of running <code>jellyfish count</code> on the Chikungunya genome",
	checks: [{
		type: "file",
		path: "chikungunya.jf",
		action: "contents",
        commandObserved: `jellyfish query chikungunya.jf ACAGTGGAC`,
		commandExpected: `echo "ACAGTGG 1"`
	}]
}];
</script>

Now it's your turn! Use Jellyfish to count all 9-mers in the Chikungunya genome. To estimate the value of the `-s` parameter, use the approach described previously.

<Exercise {criteria} />

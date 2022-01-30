<script>
import { progress } from "stores/config";
import { status } from "stores/status";
import { tutorials } from "../../../stores/tutorials";
import Listings from "./components/Listings.svelte";

const PREREQ_ID = "terminal-basics";
const PREREQ_STATUS_DONE = "done";
const PREREQ_STATUS_STARTED = "started";
const PREREQ_STATUS_UNSEEN = "unseen";
const tutorial = $tutorials.find(t => t.id === PREREQ_ID);

let prereqStatus = null;  // done, started, unseen
$: if($status.app) {
	const step = $progress[PREREQ_ID]?.step;
	if(step == null || step == 0) {
		prereqStatus = PREREQ_STATUS_UNSEEN;
	} else if(step === tutorial.steps.length - 1) {
		prereqStatus = PREREQ_STATUS_DONE;
	} else {
		prereqStatus = PREREQ_STATUS_STARTED;
	}
}
</script>

<p>Welcome to the Bash Wizardry course!</p>

<!-- Wait until app / $progress is loaded and ready to be used -->
{#if prereqStatus == null}
	Loading...
{:else}
	{#if prereqStatus === PREREQ_STATUS_DONE}
		<p>You already completed the prerequisite <strong><a target="_blank" href="/tutorials?id={PREREQ_ID}">Terminal Basics</a></strong> tutorial, so you're ready to get started with this course!</p>
	{:else}
		<p>Before we get started, it's encouraged (but not required) to finish the <strong>Terminal Basics</strong> tutorial before you proceed with this course. That tutorial will introduce you to important concepts to help you get up to speed.</p>

		{#if prereqStatus === PREREQ_STATUS_STARTED}
			<p>You already completed {$progress[PREREQ_ID]?.step} out of {tutorial.steps.length - 1} steps.</p>
		{/if}

		<Listings items={[tutorial]} colMd="12" colLg="12" colXxl="12" title={null} />
	{/if}
{/if}

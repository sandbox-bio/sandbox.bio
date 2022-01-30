<script>
import { progress } from "stores/config";
import { status } from "stores/status";
import { tutorials } from "../../../stores/tutorials";
import Listings from "./components/Listings.svelte";

const PREREQ_ID = "terminal-basics";
const PREREQ_TEXT = `Before we get started, it's encouraged (though not required) to finish the <strong>Terminal Basics</strong> tutorial before you proceed with this course. That tutorial will introduce you to important concepts to help you get up to speed.`;
const tutorial = $tutorials.find(t => t.id === PREREQ_ID);

let isLoading = true;
let isPrereqDone = false;
let isPrereqStarted = false;
let prerequisite = "";
$: if($status.app) {
	const step = $progress[PREREQ_ID]?.step;
	if(step == null || step == 0) {
		prerequisite = PREREQ_TEXT;
	} else if(step === tutorial.steps.length - 1) {
		prerequisite = `You already completed the prerequisite <strong><a target="_blank" href="/tutorials?id=${PREREQ_ID}">Terminal Basics</a></strong> tutorial, so you're ready to get started with this course!`;
	} else {
		prerequisite = PREREQ_TEXT;
		// prerequisite = `You completed ${$progress[PREREQ_ID]?.step} / ${tutorial.steps.length - 1} steps of the prerequisite <strong>Terminal Basics</strong> tutorial. ${PREREQ_TEXT}`;
	}

	// 	
	// {/if}

	isLoading = false;
}

</script>
Welcome!

<!-- Wait until app / $progress is loaded and ready to be used -->
{#if isLoading}
	Loading
{:else}
	{@html prerequisite}

	<Listings items={[tutorial]} colMd="12" colLg="12" colXxl="12" title={null} showBtn={false} />
{/if}

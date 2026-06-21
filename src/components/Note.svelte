<script>
import { Accordion, AccordionItem, Alert } from "sveltestrap";

// Svelte 3 uses 'export let' for props with default values
export let type = "practical";
export let color = "light";
export let title = "";
</script>

{#if type === "more"}
    <Accordion class="my-4 shadow">
        <AccordionItem class="border-primary">
            <div slot="header" class="fw-bold">
                {title || "🧐 Tell me more!"}
            </div>
            <div style="font-size:95%">
                <slot />
            </div>
        </AccordionItem>
    </Accordion>
{:else}
    <Alert {color} class="pb-0 border border-primary shadow my-4">
        <h6 class="fw-bold">
            {#if type === "practical"}
                🛠️ Practically speaking...
            {:else if type === "whatif"}
                🤔 What if...
            {:else if type === "interview"}
                🚀 Job interview tip
            {:else if type === "fun"}
                🙃 Fun fact
            {/if}
        </h6>
        <div style="font-size:95%">
            <slot />
        </div>
    </Alert>
{/if}

<script>
import * as marked from "marked";
let step = `Tutorial here`;

$: stepCleaned = clean(step);

// Open URLs in new tab
const renderer = {
	link(href, title, text) {
		return `<a href="${href}" target="_blank">${text}</a>`
	}
};

// Parse <Alert> tags
const tagAlert = {
	name: "alert",
	level: "block",
	start(src) { return src.match(/^<Alert>/)?.index; },
	tokenizer(src, tokens) {
		const rule = /^<Alert>[\n]*(.*)[\n]*<\/Alert>/;
		const match = rule.exec(src);
		if (match) {
			const token = {
				type: "alert",
				raw: match[1],
				text: match[1].trim(),
				tokens: []
			};
			console.log("alert", token)
			this.lexer.inline(token.text, token.tokens);
			return token;
		}
	},
	renderer(token) {
		return `
			<div class="alert alert-primary" role="alert">
				${token.text}
			</div>
		`;
	}
};

// Parse <Execute> tags
const tagExecute = {
	name: "execute",
	level: "block",
	start(src) { return src.match(/^<Execute/)?.index; },
	tokenizer(src, tokens) {
		// return false;
		const rule = /^<Execute command="(.+?)" (inline )?\/>/;
		const match = rule.exec(src);
		console.log("MATCHHHH", match)
		if (match) {
			const token = {
				type: "execute",
				raw: match[1],
				text: match[1].trim(),
				tokens: []
			};
			console.log("token", token)
			this.lexer.inline(token.text, token.tokens);
			return token;
		}
	},
	renderer(token) {
		console.log("TOKEN", token)
		// return;
		// 		<kbd on:click={exec}>
		// 		{@html command.replace(/ \\ /g, " \\ <br />&nbsp;&nbsp;&nbsp;")}
		// 	</kbd>
		return `
	<div class="card bg-dark text-white font-monospace mb-3">
		<div class="card-body">
			${token.text}
		</div>
	</div>
		`;
	}
};

marked.use({ renderer, extensions: [ tagAlert,  ] });

function clean(s) {

	return s;
}
</script>

<div class="row">
	<div class="col-6">
		<h4>Markdown</h4>
		<span class="input" role="textbox" contenteditable bind:textContent={step} />

		<!-- <textarea class="form-control" style="height: fit-content" bind:value={step}></textarea> -->
	</div>
	<div class="col-6">
		<h4>Preview</h4>
		{@html marked.parse(stepCleaned)}
	</div>
</div>

<div class="row">
	Use &lt;Alert&gt;Alert goes **here**&lt;/Alert&gt;
</div>

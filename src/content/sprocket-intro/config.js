// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	// Metadata
	id: "sprocket-intro",
	name: "Running WDL workflows with Sprocket",
	icon: "gear-wide-connected",
	date: "November 2026",
	subtitle: `by the <a href="https://sprocket.bio" target="_blank">Sprocket</a> team`,
	description: "Check, lint and run a <code>WDL</code> workflow with Sprocket, entirely in your browser.",
	tags: ["sprocket", "wdl", "workflows", "bedtools"],
	tools: ["sprocket", "bedtools"],
	difficulty: ["beginner"],

	// Lessons
	steps: [
		{ name: "WDL workflows with Sprocket", component: Intro },
		{ name: "The workflow", component: Step1 },
		{ name: "Static analysis", component: Step2, subtitle: "Check and lint", header: true },
		{ name: "Static analysis", component: Step3, subtitle: "Inputs" },
		{ name: "Running", component: Step4, subtitle: "Run the workflow", header: true },
		{ name: "Running", component: Step5, subtitle: "Explore the outputs" },
		{ name: "Running", component: Step6, subtitle: "Change the inputs" },
		{ name: "The end", component: Conclusion, header: true }
	],

	// Files within "data/" that are copied into ~/tutorial
	files: ["overlaps.wdl", "inputs.json", "sprocket.toml", "exons.bed", "cpg.bed", "fHeart-DS15839.bed", "fHeart-DS16621.bed", "fSkin-DS19745.bed"],

	// Without this, Sprocket's TLS setup reads every file in /etc/ssl/certs, each of which is
	// fetched lazily over 9p, which makes the first command take minutes
	init: "export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt"
};

// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
// import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Step7 from "./steps/Step7.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "bash-script",
	date: "July 2025",
	name: "How to write a Bash script",
	engine: "bash",
	icon: "terminal",
	subtitle: `by <a href="https://robert.bio" target="_blank">Robert Aboukhalil</a>`,
	description: "How to turn your one-liners into reusable command-line scripts",
	tags: ["terminal"],
	tools: [],
	difficulty: ["beginner"],
	tabs: [{ name: "hello.sh", contents: `echo "Hello world"\n` }],
	steps: [
		{ name: "How to write a Bash script", component: Intro, subtitle: "Why scripts?" },
		{ name: "Your first script", component: Step1, subtitle: "Hello world", header: true },
		{ name: "Your first script", component: Step2, subtitle: "Don't forget the shebang" },
		// { name: "Pardon the interruption...", component: Step3 },
		{ name: "User input", component: Step4, subtitle: "Parsing command-line arguments", header: true },
		{ name: "User input", component: Step5, subtitle: "Required parameters" },
		{ name: "User input", component: Step6, subtitle: "Optional parameters" },
		{ name: "Editing files on the command line", component: Step7, subtitle: "How to use vim", header: true },
		{ name: "The End", component: Conclusion, header: true }
	],
	files: ["hello.sh"]
};

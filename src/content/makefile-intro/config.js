// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Step7 from "./steps/Step7.md";
import Step8 from "./steps/Step8.md";
import Step9 from "./steps/Step9.md";
import Step10 from "./steps/Step10.md";
import Step11 from "./steps/Step11.md";
import Step12 from "./steps/Step12.md";
import Step13 from "./steps/Step13.md";
import Step14 from "./steps/Step14.md";
import Step15 from "./steps/Step15.md";
import Step16 from "./steps/Step16.md";
import Step17 from "./steps/Step17.md";
import Step18 from "./steps/Step18.md";
import Step19 from "./steps/Step19.md";
import Example from "./steps/Example.md";
import Conclusion from "./steps/Conclusion.md";
// Exercises
import Exercise1 from "./exercises/Exercise1.md";
import Exercise2 from "./exercises/Exercise2.md";
import Exercise3 from "./exercises/Exercise3.md";
import Exercise4 from "./exercises/Exercise4.md";

export const config = {
	// Metadata
	id: "makefile-intro",
	name: "Intro to Makefiles",
	subtitle: `by <a href="https://github.com/muppetjones" target="_blank">Stephen Bush</a>`,
	description: "Learn the basics of makefile automation.",
	tags: ["terminal", "shell", "make", "makefile", "vim"],
	difficulty: ["beginner"],

	// Preload these tools as soon as the page loads
	tools: ["make", "ls", "vim", "cat", "less", "rm"],

	// Sets up an IDE layout
	tabs: [
		{ name: "Makefile", contents: "hello-world:\n\techo \"Hello, world.\"\n" },
		{ name: "Makefile.story", contents: "story:" },
	],

	// Order of steps
	steps: [
		{ name: "An Intro to Makefiles", component: Intro },

		// Introduce makefile structure
		{ name: "Structure of a Rule", component: Step1 },

		// Recipe prefixes
		{ name: "Silence is Golden", component: Step2 },

		// PHONY vs file targets
		{ name: "Rules are made to be broken", component: Step3 },
		{ name: "Fake it 'til you make it", component: Step4 },

		// (Special) Variables
		{ name: "Bravely default", component: Step5 },

		// Prerequisites (steup)
		{ name: "Story time", component: Step6 },

		// Include directive
		{ name: "Inclusion is not an option", component: Step7 },

		// Patterns
		{ name: "Me first", component: Step8 },

		// WIldcards
		{ name: "Less is more?", component: Step9 },

		// Dynamic vs concrete rules
		{ name: "Explicitly forbidden", component: Step10 },
		{ name: "A list apart", component: Step11 },

		// Building file-based rules
		{ name: "Change is the only constant", component: Step12 },
		{ name: "Finally! File targets.", component: Step13 },
		{ name: "Aside: Make Version and Unused Variables", component: Step14 },
		{ name: "DIY functions", component: Step15 },

		// The final run
		{ name: "Almost there", component: Step16 },
		{ name: "Stay on target", component: Step17 },
		{ name: "You're all clear, kid", component: Step18 },
		{ name: "One in a million", component: Step19 },

		// Examples
		{ name: "Mad Lib", component: Example },

		{ name: "The end", component: Conclusion, header: true },

		// Exercises
		{ name: "Exercises", component: Exercise1, subtitle: "Gzip text files", header: true },
		{ name: "Exercises", component: Exercise2, subtitle: "Ungzip text files" },
		{ name: "Exercises", component: Exercise3, subtitle: "Convert and index SAM files" },
		{ name: "Exercises", component: Exercise4, subtitle: "Madlib!" },
	],


	// Files within "data/" that you need at runtime.
	// Whenever you update files within "data/", you will need to restart the "./setup.sh" script.
	files: [
		"Makefile.madlib",
		"generate_word.sh",
		"generate_madlib_part.sh"
	]
};

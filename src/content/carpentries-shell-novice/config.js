// Steps
import Intro from "./steps/Intro.md";
import Step2_Exercise1 from "./steps/Step2_Exercise1.md";
import Step2_Exercise2 from "./steps/Step2_Exercise2.md";
import Step2_Exercise3 from "./steps/Step2_Exercise3.md";
import Step2_Summary from "./steps/Step2_Summary.md";
import Step3_Summary from "./steps/Step3_Summary.md";
import Step3_Exercise1 from "./steps/Step3_Exercise1.md";
import Step3_Exercise2 from "./steps/Step3_Exercise2.md";
import Step3_Exercise3 from "./steps/Step3_Exercise3.md";
import Step3_Exercise4 from "./steps/Step3_Exercise4.md";
import Step4_Exercise1 from "./steps/Step4_Exercise1.md";
import Step4_Exercise2 from "./steps/Step4_Exercise2.md";
import Step4_Exercise3 from "./steps/Step4_Exercise3.md";
import Step4_Summary from "./steps/Step4_Summary.md";
import Step5_Summary from "./steps/Step5_Summary.md";
import Step6_Summary from "./steps/Step6_Summary.md";
import Step7_Summary from "./steps/Step7_Summary.md";

export const config = {
	id: "carpentries-shell-novice",
	listed: false,
	name: "The Unix Shell",
	subtitle: `by <a href="https://swcarpentry.github.io/shell-novice/04-pipefilter.html" target="_blank">The Software Carpentry</a>`,
	description: "Example Carpentries tutorial",
	tags: ["unix", "shell", "terminal"],
	difficulty: ["beginner"],
	init: "mkdir -p ~/tutorial/raw",
	steps: [
		{ name: "The Unix Shell", component: Intro },
		{ name: "Navigating Files and Directories", component: Step2_Exercise1, subtitle: "Exercise 1: Absolute vs Relative Paths", header: true },
		{ name: "Navigating Files and Directories", component: Step2_Exercise2, subtitle: "Exercise 2: Absolute vs Relative Paths" },
		{ name: "Navigating Files and Directories", component: Step2_Exercise3, subtitle: "Exercise 3: ls Reading Comprehension" },
		{ name: "Navigating Files and Directories", component: Step2_Summary, subtitle: "Key Points" },
		{ name: "Working With Files and Directories", component: Step3_Exercise1, subtitle: "Exercise 1: Moving Files to a new folder", header: true },
		{ name: "Working With Files and Directories", component: Step3_Exercise2, subtitle: "Exercise 2: Renaming Files" },
		{ name: "Working With Files and Directories", component: Step3_Exercise3, subtitle: "Exercise 3: Moving and Copying" },
		{ name: "Working With Files and Directories", component: Step3_Exercise4, subtitle: "Exercise 4: List filenames matching a pattern" },
		{ name: "Working With Files and Directories", component: Step3_Summary, subtitle: "Key Points" },
		{ name: "Pipes and Filters", component: Step4_Exercise1, subtitle: "Exercise 1: Appending Data", header: true },
		{ name: "Pipes and Filters", component: Step4_Exercise2, subtitle: "Exercise 2: Piping Commands Together" },
		{ name: "Pipes and Filters", component: Step4_Exercise3, subtitle: "Exercise 3: Pipe Construction" },
		{ name: "Pipes and Filters", component: Step4_Summary, subtitle: "Key Points" },
		{ name: "Loops", component: Step5_Summary, subtitle: "Key Points", header: true },
		{ name: "Shell Scripts", component: Step6_Summary, subtitle: "Key Points", header: true },
		{ name: "Finding Things", component: Step7_Summary, subtitle: "Key Points", header: true }
	],
	files: [
		"analyzed/fructose.dat",
		"analyzed/glucose.dat",
		"analyzed/maltose.dat",
		"analyzed/sucrose.dat",
		"exercise-data/numbers.txt",
		"exercise-data/alkanes/propane.pdb",
		"exercise-data/alkanes/octane.pdb",
		"exercise-data/alkanes/cubane.pdb",
		"exercise-data/alkanes/ethane.pdb",
		"exercise-data/alkanes/pentane.pdb",
		"exercise-data/alkanes/methane.pdb",
		"exercise-data/animal-counts/animals.csv",
		"exercise-data/creatures/minotaur.dat",
		"exercise-data/creatures/unicorn.dat",
		"exercise-data/creatures/basilisk.dat",
		"exercise-data/writing/LittleWomen.txt",
		"exercise-data/writing/haiku.txt",
		"north-pacific-gyre/NENE01729B.txt",
		"north-pacific-gyre/NENE02040B.txt",
		"north-pacific-gyre/NENE01729A.txt",
		"north-pacific-gyre/NENE02040A.txt",
		"north-pacific-gyre/NENE01971Z.txt",
		"north-pacific-gyre/goodiff.sh",
		"north-pacific-gyre/NENE01978B.txt",
		"north-pacific-gyre/NENE01978A.txt",
		"north-pacific-gyre/NENE01812A.txt",
		"north-pacific-gyre/NENE01736A.txt",
		"north-pacific-gyre/NENE01843B.txt",
		"north-pacific-gyre/goostats.sh",
		"north-pacific-gyre/NENE01843A.txt",
		"north-pacific-gyre/NENE02018B.txt",
		"north-pacific-gyre/NENE02040Z.txt",
		"north-pacific-gyre/NENE02043A.txt",
		"north-pacific-gyre/NENE01751A.txt",
		"north-pacific-gyre/NENE02043B.txt",
		"north-pacific-gyre/NENE01751B.txt"
	]
};

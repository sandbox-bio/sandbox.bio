// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";

export const config = {
	id: "carpentries-shell-novice",
	listed: false,
	name: "The Unix Shell",
	subtitle: `by <a href="https://swcarpentry.github.io/shell-novice/04-pipefilter.html" target="_blank">The Software Carpentry</a>`,
	description: "Example Carpentries tutorial",
	tags: ["unix", "shell", "terminal"],
	difficulty: ["beginner"],
	steps: [
		{ name: "The Unix Shell", component: Intro },
		{ name: "Pipes and Filters", component: Step1 },
		{ name: "Loops", component: Step2 },
		{ name: "Shell Scripts", component: Step3 }
	],
	files: [
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

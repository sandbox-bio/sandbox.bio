// Steps
import Intro from "./steps/Intro.md";
import Episode6_Exercise1 from "./steps/Episode6_Exercise1.md";
import Episode6_Exercise2 from "./steps/Episode6_Exercise2.md";
import Episode6_Exercise3 from "./steps/Episode6_Exercise3.md";
import Episode7_Exercise1 from "./steps/Episode7_Exercise1.md";
import Episode7_Exercise2 from "./steps/Episode7_Exercise2.md";
import Episode7_Exercise3 from "./steps/Episode7_Exercise3.md";
import Episode7_Exercise4 from "./steps/Episode7_Exercise4.md";

export const config = {
	id: "carpentries-shell-novice",
	listed: false,
	name: "Terminal Exercises",
	subtitle: "Exercises from the Carpentries' Unix Shell lesson",
	description: "Exercises from the Carpentries' Unix Shell lesson",
	tags: ["terminal", "exercises", "carpentries"],
	tools: ["cut", "sort", "uniq", "grep", "cut", "wc", "vim", "find"],
	difficulty: ["beginner"],
	init: "mkdir -p ~/tutorial/raw",
	steps: [
		{ name: "Terminal Exercises", component: Intro },
		{ name: "Shell Scripts", component: Episode6_Exercise1, subtitle: "List unique species", header: true },
		{ name: "Shell Scripts", component: Episode6_Exercise2, subtitle: "Variables in shell scripts" },
		{ name: "Shell Scripts", component: Episode6_Exercise3, subtitle: "Longest file with a given extension" },
		{ name: "Finding Things", component: Episode7_Exercise1, subtitle: "Using grep", header: true },
		{ name: "Finding Things", component: Episode7_Exercise2, subtitle: "Tracking a species" },
		{ name: "Finding Things", component: Episode7_Exercise3, subtitle: "Little women" },
		{ name: "Finding Things", component: Episode7_Exercise4, subtitle: "Matching and subtracting" }
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

// Steps
import Intro from "./steps/Intro.md";
import Episode2_Exercise1 from "./steps/Episode2_Exercise1.md";
import Episode2_Exercise2 from "./steps/Episode2_Exercise2.md";
import Episode2_Exercise3 from "./steps/Episode2_Exercise3.md";
import Episode3_Exercise1 from "./steps/Episode3_Exercise1.md";
import Episode3_Exercise2 from "./steps/Episode3_Exercise2.md";
import Episode3_Exercise3 from "./steps/Episode3_Exercise3.md";
import Episode3_Exercise4 from "./steps/Episode3_Exercise4.md";
import Episode4_Exercise1 from "./steps/Episode4_Exercise1.md";
import Episode4_Exercise2 from "./steps/Episode4_Exercise2.md";
import Episode4_Exercise3 from "./steps/Episode4_Exercise3.md";
import Episode4_Exercise4 from "./steps/Episode4_Exercise4.md";
import Episode4_Exercise5 from "./steps/Episode4_Exercise5.md";
import Episode5_Exercise1 from "./steps/Episode5_Exercise1.md";
import Episode5_Exercise2 from "./steps/Episode5_Exercise2.md";
import Episode5_Exercise3 from "./steps/Episode5_Exercise3.md";
import Episode5_Exercise4 from "./steps/Episode5_Exercise4.md";
import Episode6_Exercise1 from "./steps/Episode6_Exercise1.md";
import Episode6_Exercise2 from "./steps/Episode6_Exercise2.md";
import Episode6_Exercise3 from "./steps/Episode6_Exercise3.md";
import Episode7_Exercise1 from "./steps/Episode7_Exercise1.md";
import Episode7_Exercise2 from "./steps/Episode7_Exercise2.md";
import Episode7_Exercise3 from "./steps/Episode7_Exercise3.md";
import Episode7_Exercise4 from "./steps/Episode7_Exercise4.md";

export const config = {
	id: "carpentries-shell-novice",
	name: "Terminal Exercises",
	icon: "terminal-plus",
	subtitle: "Exercises from the Carpentries' Unix Shell lesson",
	description: "Command-line exercises from the Carpentries' Unix Shell lesson.",
	tags: ["terminal", "exercises", "carpentries"],
	tools: ["cut", "sort", "uniq", "grep", "cut", "wc", "vim", "find"],
	difficulty: ["intermediate"],
	new: true,
	init: "mkdir -p $TUTORIAL/raw $TUTORIAL/loops",
	steps: [
		{ name: "Terminal Exercises", component: Intro },
		// ---------------------------------------------------------------------
		// Episode 2
		// ---------------------------------------------------------------------
		// https://swcarpentry.github.io/shell-novice/02-filedir.html#exploring-more-ls-options
		// -
		// https://swcarpentry.github.io/shell-novice/02-filedir.html#listing-in-reverse-chronological-order
		// -
		// https://swcarpentry.github.io/shell-novice/02-filedir.html#absolute-vs-relative-paths
		{ name: "Navigating Files and Directories", component: Episode2_Exercise1, subtitle: "Absolute vs relative paths", header: true },
		// https://swcarpentry.github.io/shell-novice/02-filedir.html#relative-path-resolution
		{ name: "Navigating Files and Directories", component: Episode2_Exercise2, subtitle: "Absolute vs relative paths" },
		// https://swcarpentry.github.io/shell-novice/02-filedir.html#ls-reading-comprehension
		{ name: "Navigating Files and Directories", component: Episode2_Exercise3, subtitle: "Ls reading comprehension" },

		// ---------------------------------------------------------------------
		// Episode 3
		// ---------------------------------------------------------------------
		// https://swcarpentry.github.io/shell-novice/03-create.html#creating-files-a-different-way
		// -
		// https://swcarpentry.github.io/shell-novice/03-create.html#moving-files-to-a-new-folder
		{ name: "Working With Files and Directories", component: Episode3_Exercise1, subtitle: "Moving Files to a new folder", header: true },
		// https://swcarpentry.github.io/shell-novice/03-create.html#renaming-files
		{ name: "Working With Files and Directories", component: Episode3_Exercise2, subtitle: "Renaming Files" },
		// https://swcarpentry.github.io/shell-novice/03-create.html#moving-and-copying
		{ name: "Working With Files and Directories", component: Episode3_Exercise3, subtitle: "Moving and Copying" },
		// https://swcarpentry.github.io/shell-novice/03-create.html#using-rm-safely
		// -
		// https://swcarpentry.github.io/shell-novice/03-create.html#copy-with-multiple-filenames
		// -
		// https://swcarpentry.github.io/shell-novice/03-create.html#list-filenames-matching-a-pattern
		{ name: "Working With Files and Directories", component: Episode3_Exercise4, subtitle: "List filenames matching a pattern" },
		// https://swcarpentry.github.io/shell-novice/03-create.html#more-on-wildcards
		// -
		// https://swcarpentry.github.io/shell-novice/03-create.html#organizing-directories-and-files
		// -
		// https://swcarpentry.github.io/shell-novice/03-create.html#reproduce-a-folder-structure
		// -

		// ---------------------------------------------------------------------
		// Episode 4
		// ---------------------------------------------------------------------
		// https://swcarpentry.github.io/shell-novice/04-pipefilter.html#what-does-sort--n-do
		// -
		// https://swcarpentry.github.io/shell-novice/04-pipefilter.html#what-does-mean
		// -
		// https://swcarpentry.github.io/shell-novice/04-pipefilter.html#appending-data
		{ name: "Pipes and Filters", component: Episode4_Exercise1, subtitle: "Appending data", header: true },
		// https://swcarpentry.github.io/shell-novice/04-pipefilter.html#piping-commands-together
		{ name: "Pipes and Filters", component: Episode4_Exercise2, subtitle: "Piping commands together" },
		// https://swcarpentry.github.io/shell-novice/04-pipefilter.html#pipe-reading-comprehension
		// -
		// https://swcarpentry.github.io/shell-novice/04-pipefilter.html#pipe-construction
		{ name: "Pipes and Filters", component: Episode4_Exercise3, subtitle: "Pipe construction" },
		// https://swcarpentry.github.io/shell-novice/04-pipefilter.html#which-pipe
		{ name: "Pipes and Filters", component: Episode4_Exercise4, subtitle: "Which pipe?" },
		// https://swcarpentry.github.io/shell-novice/04-pipefilter.html#removing-unneeded-files
		{ name: "Pipes and Filters", component: Episode4_Exercise5, subtitle: "Removing unneeded files" },

		// ---------------------------------------------------------------------
		// Episode 5
		// ---------------------------------------------------------------------
		// https://swcarpentry.github.io/shell-novice/05-loop.html#write-your-own-loop
		{ name: "Loops", component: Episode5_Exercise1, subtitle: "Write your own loop", header: true },
		// https://swcarpentry.github.io/shell-novice/05-loop.html#variables-in-loops
		// -
		// https://swcarpentry.github.io/shell-novice/05-loop.html#limiting-sets-of-files
		// https://swcarpentry.github.io/shell-novice/05-loop.html#limiting-sets-of-files CONTINUED
		{ name: "Loops", component: Episode5_Exercise2, subtitle: "Limiting sets of files" },
		// https://swcarpentry.github.io/shell-novice/05-loop.html#saving-to-a-file-in-a-loop---part-one
		{ name: "Loops", component: Episode5_Exercise3, subtitle: "Saving to a file in a loop" },
		// https://swcarpentry.github.io/shell-novice/05-loop.html#saving-to-a-file-in-a-loop---part-two
		{ name: "Loops", component: Episode5_Exercise4, subtitle: "Saving to a file in a loop - part 2" },
		// https://swcarpentry.github.io/shell-novice/05-loop.html#doing-a-dry-run
		// -
		// https://swcarpentry.github.io/shell-novice/05-loop.html#nested-loops
		// -

		// ---------------------------------------------------------------------
		// Episode 6
		// ---------------------------------------------------------------------
		// https://swcarpentry.github.io/shell-novice/06-script.html#list-unique-species
		{ name: "Shell Scripts", component: Episode6_Exercise1, subtitle: "List unique species", header: true },
		// https://swcarpentry.github.io/shell-novice/06-script.html#why-record-commands-in-the-history-before-running-them
		// -
		// https://swcarpentry.github.io/shell-novice/06-script.html#variables-in-shell-scripts
		{ name: "Shell Scripts", component: Episode6_Exercise2, subtitle: "Variables in shell scripts" },
		// https://swcarpentry.github.io/shell-novice/06-script.html#find-the-longest-file-with-a-given-extension
		{ name: "Shell Scripts", component: Episode6_Exercise3, subtitle: "Longest file with a given extension" },
		// https://swcarpentry.github.io/shell-novice/06-script.html#script-reading-comprehension
		// -
		// https://swcarpentry.github.io/shell-novice/06-script.html#debugging-scripts
		// -

		// ---------------------------------------------------------------------
		// Episode 7
		// ---------------------------------------------------------------------
		// https://swcarpentry.github.io/shell-novice/07-find.html#using-grep
		{ name: "Finding Things", component: Episode7_Exercise1, subtitle: "Using grep", header: true },
		// https://swcarpentry.github.io/shell-novice/07-find.html#tracking-a-species
		{ name: "Finding Things", component: Episode7_Exercise2, subtitle: "Tracking a species" },
		// https://swcarpentry.github.io/shell-novice/07-find.html#little-women
		{ name: "Finding Things", component: Episode7_Exercise3, subtitle: "Little women" },
		// https://swcarpentry.github.io/shell-novice/07-find.html#matching-and-subtracting
		{ name: "Finding Things", component: Episode7_Exercise4, subtitle: "Matching and subtracting" }
		// https://swcarpentry.github.io/shell-novice/07-find.html#find-pipeline-reading-comprehension
		// -
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

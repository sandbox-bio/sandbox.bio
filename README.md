# sandbox.bio

Interactive bioinformatics command-line tutorials for bioinformatics at [sandbox.bio](https://sandbox.bio).

# Contributing tutorials

## Overview

To contribute a new tutorial:

1. Reach out through [GitHub Discussions](https://github.com/sandbox-bio/sandbox.bio/discussions) or by email to discuss the tutorial you want to add to sandbox.bio. It can be a new tutorial you want to write or an existing tutorial you want to port over. Please include the name of the command-line tools your tutorial needs.
2. We'll work together to add new tools to the sandbox.bio engine's [Dockerfile](https://github.com/sandbox-bio/v86/blob/master/tools/docker/debian/Dockerfile). To be supported, tools must be compiled to a 32-bit i686 architecture, and if they use SIMD, can only use SSE, SSE2, and SSE3. JVM-based tools are not supported because they require large downloads, and are not performant in this environment.
3. Fork this repo and set up your local environment (next section)
4. Write the tutorial
5. Send a PR

## Local development setup

Before setting up your local environment, make sure you:

- Have forked this repo so you can send a PR to contribute your changes.
- Have `npm` and `node` [installed on your computer](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm).

Then, on your terminal, run the following commands:

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME_GOES_HERE/sandbox.bio.git
cd sandbox.bio
git checkout -b my-new-tutorial

# Set up your environment
./bin/setup.sh
```

If everything goes well, a browser window will open with sandbox.bio running locally, with a blue "Development Mode" ribbon at the top of the page. You may see a blank page for a few seconds.

## Writing your tutorial

Make a copy of the template:

```bash
TUTORIAL_ID=awesome-tool-intro
mkdir -p "./src/content/$TUTORIAL_ID"
cp -r "./src/content/_template/" "./src/content/$TUTORIAL_ID"
```

Then you need to modify `src/stores/tutorials.js`:

1. Import your new tutorial at the [top of the file](https://github.com/sandbox-bio/sandbox.bio/blob/main/src/stores/tutorials.js#L30): `import { config as awesomeToolIntro } from "$content/awesome-tool-intro/config";`
2. Add your tutorial to the `tutorials` variable on [this line](https://github.com/sandbox-bio/sandbox.bio/blob/main/src/stores/tutorials.js#L54)

Modify the files as needed--see the next section for information about how files are organized.

### Tutorial Structure

Here's how you should structure your tutorial within the sandbox.bio repository:

```bash
/src/content/
  some-tool-intro/
    data/               # Data you want to be preloaded in the terminal when the tutorial loads.
      reads.fastq       # Keep these files as small as possible, ideally < 100KB if possible.
      aligned.bam
    steps/              # Tutorial content for each step, defined in Markdown format (see format below).
        Step1.md
        Step2.md
    exercises/          # Exercises (see format below).
      Exercise1.md
    config.js           # Configuration file that defines order of tutorial steps, and other metadata.
    README.md           # Include where you downloaded files from, and how/if they were processed (optional).
```

### Tutorial format

Tutorial steps are defined as Markdown. See the [bedtools tutorial](https://raw.githubusercontent.com/sandbox-bio/sandbox.bio/main/src/content/bedtools-intro/steps/Step12.md) for an example format.

### Exercises

See [this exercise](https://raw.githubusercontent.com/sandbox-bio/sandbox.bio/b3174e01e25c48c1bf655e894626eb0a09c88992/src/content/debugging-puzzles/steps/PuzzleBedSpaces.md) for how to format an exercise, optional hints to show learners, and the criteria for marking an exercise as successfully completed.

### Quizzes

Multiple-choice quizzes can be implemented as follows. If there is only one valid answer, the quiz will display radio boxes, otherwise it will show checkboxes.

```js
<Quiz
	id="step3-quiz2"
	choices={[
		{ valid: false, value: `Montreal` },
		{ valid: true, value: `Ottawa` },
		{ valid: false, value: `Toronto` },
		{ valid: false, value: `Vancouver` }
	]}
>
	<span slot="prompt">What is the capital of Canada?</span>
</Quiz>
```

### config.js

See the [bedtools tutorial](https://github.com/sandbox-bio/sandbox.bio/blob/main/src/content/bedtools-intro/config.js) for an example `config.js` to use as baseline.

---

# Gists (beta)

You can use the sandbox.bio command line environment to share small snippets of code, make them interactive, and include sample data files:

|GitHub Gist||sandbox.bio|
|--|--|--|
| <img src="https://github.com/user-attachments/assets/81129198-aeba-417f-a080-bca31153d93c" height="320"> | ➡️ | <img src="https://github.com/user-attachments/assets/a556c66e-46ba-4760-9c11-8b8ceaca06c7" height="400"> |
|[gist.github.com/9ee9bd2dd6be673acdc6971e852a28a3](https://gist.github.com/robertaboukhalil/9ee9bd2dd6be673acdc6971e852a28a3)||[sandbox.bio/gists/9ee9bd2dd6be673acdc6971e852a28a3](https://sandbox.bio/gists/9ee9bd2dd6be673acdc6971e852a28a3)|

1. Create a GitHub Gist at https://gist.github.com/
	- [ ] Create a `.md` Markdown file with notes about the commands. Bash code blocks are runnable within sandbox.bio
	- [ ] Add data files if relevant. These **files will be available in the sandbox.bio command-line environment at load time**!
2. Once created, the URL contains the gist ID, e.g. gist.github.com/gists/username/gist_id
3. All you have to do is go to sandbox.bio/gists/gist_id (no need for the username, just the gist ID)

sandbox.bio will render the Markdown file(s) from the Gist, show code blocks with syntax highlighting, and allow users to run the code blocks in their sandbox.bio terminal environment, where data files will be ready for them.

Limitations:
* GitHub Gists don't support binary files except common ones like zip files
* For larger files, you can use the GitHub command line tool to upload them more easily:
  * Upload local file: `gh gist edit GIST_ID_HERE ./path/to/data.fastq --add data.fastq`
  * Modify existing file: `gh gist edit GIST_ID_HERE ./path/to/data.fastq --filename data.fastq`

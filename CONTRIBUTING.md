# Contributing Tutorials

### Overview

Here are the steps to contribute a new tutorial:

1. Reach out through [GitHub Discussions](https://github.com/sandbox-bio/sandbox.bio/discussions) or email to discuss the tutorial you want to add to sandbox.bio. It can be a new tutorial you want to write or an existing tutorial you want to port over.
2. If needed, we'll work together to add new tools to the sandbox.bio engine's [Dockerfile](https://github.com/sandbox-bio/v86/blob/master/tools/docker/debian/Dockerfile). To be supported, tools must be compiled to a 32-bit i686 architecture, and if they use SIMD, can only use SSE, SSE2, and SSE3. JVM-based tools are not supported because they require large downloads, and are not performant in this environment.
3. Draft the tutorial by creating a Markdown file for each step of the tutorial (details below).
4. Fork the repo, and send a Pull Request to merge your branch into `main`.

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

### Steps

Tutorial steps are defined as Markdown. See the [bedtools tutorial](https://raw.githubusercontent.com/sandbox-bio/sandbox.bio/main/src/content/bedtools-intro/steps/Step12.md) for an example format.

Make sure to check out [sandbox.bio Studio](https://sandbox.bio/studio) for help designing tutorial steps.

### Exercises

See [this exercise](https://raw.githubusercontent.com/sandbox-bio/sandbox.bio/b3174e01e25c48c1bf655e894626eb0a09c88992/src/content/debugging-puzzles/steps/PuzzleBedSpaces.md) for how to format an exercise, optional hints to show learners, and the criteria for marking an exercise as successfully completed.

### Quizzes

Multiple-choice quizzes can be implemented as follows. If there is only one valid answer, the quiz will display radio boxes, otherwise it will show checkboxes.

```js
<Quiz id="step3-quiz2" choices={[
    { valid: false, value: `Montreal`},
    { valid: true, value: `Ottawa`},
    { valid: false, value: `Toronto`},
    { valid: false, value: `Vancouver`},
]}>
    <span slot="prompt">
        What is the capital of Canada?
    </span>
</Quiz>
```

### config.js

See the [bedtools tutorial](https://github.com/sandbox-bio/sandbox.bio/blob/main/src/content/bedtools-intro/config.js) for an example `config.js` to use as baseline.

# sandbox.bio

Interactive bioinformatics command-line tutorials.

## Contributing

### Overview

To contribute a new tutorial:

1. Reach out through [GitHub Discussions](https://github.com/sandbox-bio/sandbox.bio/discussions) or by email to discuss the tutorial you want to add to sandbox.bio. It can be a new tutorial you want to write or an existing tutorial you want to port over. Please include the name of the command-line tools your tutorial needs.
2. We'll work together to add new tools to the sandbox.bio engine's [Dockerfile](https://github.com/sandbox-bio/v86/blob/master/tools/docker/debian/Dockerfile). To be supported, tools must be compiled to a 32-bit i686 architecture, and if they use SIMD, can only use SSE, SSE2, and SSE3. JVM-based tools are not supported because they require large downloads, and are not performant in this environment.
3. Set up your local environment (next section)
4. Write the tutorial

### Local development setup

Before setting up your local environment, make sure you have `npm` and `node` [installed on your computer](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm).

Then, run the following commands:

```bash
git clone https://github.com/sandbox-bio/sandbox.bio.git
cd sandbox.bio
./bin/setup.sh
```

If everything goes well, a browser window will open with sandbox.bio running locally, with a blue "Development Mode" ribbon at the top of the page. You may see a blank page for a few seconds.


# sandbox.bio

Interactive bioinformatics command-line tutorials.

## Contributing

### Overview

To contribute a new tutorial:

1. Reach out through [GitHub Discussions](https://github.com/sandbox-bio/sandbox.bio/discussions) or email to discuss the tutorial you want to add to sandbox.bio. It can be a new tutorial you want to write or an existing tutorial you want to port over. Please include the name of the tools your tutorial will need.
2. If needed, we'll work together to add new tools to the sandbox.bio engine's [Dockerfile](https://github.com/sandbox-bio/v86/blob/master/tools/docker/debian/Dockerfile). To be supported, tools must be compiled to a 32-bit i686 architecture, and if they use SIMD, can only use SSE, SSE2, and SSE3. JVM-based tools are not supported because they require large downloads, and are not performant in this environment.
3. Write your tutorial (see next section)


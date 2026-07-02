<script>
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

In this tutorial, we'll cover the basics of using a Makefile. If you're familiar with
Snakemake, then you're already familiar with the basic concept of a make file. Although
`make` predates Snakemake by nearly 40 years, the basic premise is the same: Define a
recipe under a target, and execute the recipe by calling the target from the command
line.

Originally, `make` was developed to streamline the compilation of C code (don't worry,
we're not going to force you into C), and many newer tools aim to provide the same
functionality in a newer package. However, `make` remains a simple, stable staple across
many software projects.

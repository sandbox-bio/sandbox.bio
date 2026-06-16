<script>
import Link from "$components/Link.svelte";
</script>

This concludes the `make` tutorial!

You now know enough to manage your own C++ project! (that was your goal, right?)

Even if you don't ever plan on programming in C or another compiled language, `make` is
a powerful tool in your programming belt. It is used across many open source projects,
and it behaves the same regardless of the language of the project you're using. With the
Windows Subsystem for Linux (WSL), `make` is also OS-agnostic.

Here's a couple of practical exercises you can try on your own:

- Write recipes to compress or decompress any FASTQ files in your project.
- Write recipes to convert SAM files to BAM -- and index them.

Here are additional resources to check out:

- <Link href="https://man7.org/linux/man-pages/man1/make.1.html">Official manpage for `make`</Link>
- <Link href="https://www.gnu.org/software/make/manual/html_node/index.html">Official `make` docs</Link> (may be broken)
- <Link href="https://makefiletutorial.com/">Makefile tutorial</Link>
- <Link href="https://github.com/arq5x/bedtools2">Codebase</Link>

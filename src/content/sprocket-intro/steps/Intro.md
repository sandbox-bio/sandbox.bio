<script>
import Link from "$components/Link.svelte";
</script>

<Link href="https://sprocket.bio">Sprocket</Link> is a command-line tool for working with the <Link href="https://openwdl.org">Workflow Description Language</Link> (WDL). It checks, lints, formats and runs WDL workflows.

Everything in this tutorial runs **inside your browser**. The terminal on the right is a small Linux machine emulated in WebAssembly, with `sprocket` and `bedtools` already installed. You don't need to install anything.

We'll write nothing from scratch. Instead, we'll take a small workflow that counts how many DNase I hypersensitivity peaks overlap exons, using data from the <Link href="/tutorials/bedtools-intro">bedtools tutorial</Link>, and:

1. Read the workflow
2. Check and lint it
3. Validate its inputs
4. Run it and explore the results

> The sandbox has no Docker, so Sprocket runs each task directly on the machine. On your laptop or cluster, the same WDL runs in containers without changes.

<script>
import Execute from "$components/Execute.svelte";
</script>

Alongside the text report, `fastp` also creates a human-readable HTML report (`fastp.html`), and a JSON report (`fastp.json`) meant to be analyzed by a program (more on that in a bit): <Execute inline command="ls -l" />

Let's open up this HTML report in a new tab:

<Execute command="open fastp.html" />

At the top of the HTML file, you'll find the same stats we saw in the previous step, but the rest of the file contains interactive plots showing, for example, how the base quality changes as a function of the position along the read.

Note that the y-axis does not start at 0 so the decrease in quality looks more pronounced than it actually is; in fact, the base quality is > 30 across the length of the read, which is excellent!

> The `open` command will also work on MacOS, but on Linux systems, you can use `xdg-open` instead.

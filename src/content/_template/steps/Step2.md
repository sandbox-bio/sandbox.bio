<script>
import Execute from "$components/Execute.svelte";
</script>

Bedtools is a command-line tool. To bring up the help, just type <Execute command={"bedtools"} inline />

As you can see, there are multiple "subcommands" and for `bedtools` to work you must tell it which subcommand you want to use.

Examples:

- <Execute command={"bedtools intersect"} inline />
- <Execute command={"bedtools merge"} inline />
- <Execute command={"bedtools subtract"} inline />

What version am I using?:

- <Execute command={"bedtools --version"} inline />

How can I get more help?:

- <Execute command={"bedtools --contact"} inline />

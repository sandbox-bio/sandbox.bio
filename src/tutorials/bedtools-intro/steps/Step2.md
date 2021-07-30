<script>
import Execute from "../../Execute.svelte";
</script>

Bedtools is a command-line tool. To bring up the help, just type <Execute command={"bedtools"} />

As you can see, there are multiple "subcommands" and for `bedtools` to work you must tell it which subcommand you want to use.

Examples: 

* <Execute command={"bedtools intersect"} />
* <Execute command={"bedtools merge"} />
* <Execute command={"bedtools subtract"} />

What version am I using?: 

* <Execute command={"bedtools --version"} />

How can I get more help?: 

* <Execute command={"bedtools --contact"} />

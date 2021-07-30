<script>
import Execute from "../../Execute.svelte";
</script>

Bedtools is a command-line tool. To bring up the help, just type <Execute command={"bedtools"} />

As you can see, there are multiple "subcommands" and for `bedtools` to work you must tell it which subcommand you want to use.

<span>
	Examples: 
	<ul>
		<li><Execute command={"bedtools intersect"} /></li>
		<li><Execute command={"bedtools merge"} /></li>
		<li><Execute command={"bedtools subtract"} /></li>
	</ul>
</span>

<span>
	What version am I using?: 
	<ul>
		<li><Execute command={"bedtools --version"} /></li>
	</ul>
</span>

<span>
	How can I get more help?: 
	<ul>
		<li><Execute command={"bedtools --contact"} /></li>
	</ul>
</span>

<script>
import Execute from "../../Execute.svelte";
</script>

<p>Bedtools is a command-line tool. To bring up the help, just type <Execute command={"bedtools"} /></p>

<p>As you can see, there are multiple "subcommands" and for <code>bedtools</code> to work you must tell it which subcommand you want to use.</p>

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

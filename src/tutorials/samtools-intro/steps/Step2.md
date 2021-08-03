<script>
import Execute from "components/Execute.svelte";
</script>

To bring up the help, just type <Execute command={"samtools"} inline={true} />

As you can see, there are multiple "subcommands" and for samtools to work, you must tell it which subcommand you want to use.

Examples: 

* <Execute command={"samtools view"} inline={true} />
* <Execute command={"samtools sort"} inline={true} />
* <Execute command={"samtools depth"} inline={true} />

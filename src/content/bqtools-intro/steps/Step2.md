<script>
import Execute from "$components/Execute.svelte";
</script>

`bqtools` is a command-line tool.
To bring up the help screen, type <Execute command="bqtools --help" inline />

As you can see, there are multiple "subcommands" and for `bqtools` to work you must tell it which subcommand you want to use.

Examples:

- <Execute command="bqtools encode --help" inline />
- <Execute command="bqtools decode --help" inline />
- <Execute command="bqtools grep --help" inline />

What version am I using?:

- <Execute command="bqtools --version" inline />

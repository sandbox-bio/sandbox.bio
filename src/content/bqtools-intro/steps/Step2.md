<script>
import Execute from "$components/Execute.svelte";
</script>

`bqtools` is a command-line tool.
To bring up the help use type `<Execute command="bqtools --help" inline />`

As you can see, there are multiple "subcommands" and for `bqtools` to work you must tell it which subcommand you want to use.

Examples:

- <Execute command="bqtools encode" inline />
- <Execute command="bqtools decode" inline />
- <Execute command="bqtools grep" inline />

What version am I using?:

- <Execute command="bqtools --version" inline />

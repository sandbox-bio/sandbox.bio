<script>
import Alert from "$components/Alert.svelte";
import Link from "$components/Link.svelte";
</script>

<Alert>
    - Save commands in files (usually called shell scripts) for re-use.
    - `bash [filename]` runs the commands saved in a file.
    - `$@` refers to all of a shell script's command-line arguments.
    - `$1`, `$2`, etc., refer to the first command-line argument, the second command-line argument, etc.
    - Place variables in quotes if the values might have spaces in them.
    - Letting users decide what files to process is more flexible and more consistent with built-in Unix commands.

    <Link href="https://swcarpentry.github.io/shell-novice/06-script.html">Visit original Carpentries tutorial</Link>

</Alert>

<script>
import Alert from "$components/Alert.svelte";
import Link from "$components/Link.svelte";
</script>

<Alert>
    - `cp [old] [new]` copies a file.
    - `mkdir [path]` creates a new directory.
    - `mv [old] [new]` moves (renames) a file or directory.
    - `rm [path]` removes (deletes) a file.
    - `*` matches zero or more characters in a filename, so `*.txt` matches all files ending in `.txt`.
    - `?` matches any single character in a filename, so `?.txt` matches `a.txt` but not `any.txt`.
    - Use of the Control key may be described in many ways, including `Ctrl-X`, `Control-X`, and `^X`.
    - The shell does not have a trash bin: once something is deleted, it's really gone.
    - Most files' names are `something.extension`. The extension isn't required, and doesn't guarantee anything, but is normally used to indicate the type of data in the file.
    - Depending on the type of work you do, you may need a more powerful text editor than Nano.

    <Link href="https://swcarpentry.github.io/shell-novice/03-create.html">Visit original Carpentries tutorial</Link>

</Alert>

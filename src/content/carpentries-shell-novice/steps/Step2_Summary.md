<script>
import Alert from "$components/Alert.svelte";
import Link from "$components/Link.svelte";
</script>

<Alert>
    - The file system is responsible for managing information on the disk.
    - Information is stored in files, which are stored in directories (folders).
    - Directories can also store other directories, which then form a directory tree.
    - `pwd` prints the user's current working directory.
    - `ls [path]` prints a listing of a specific file or directory; `ls` on its own lists the current working directory.
    - `cd [path]` changes the current working directory.
    - Most commands take options that begin with a single `-`.
    - Directory names in a path are separated with `/` on Unix, but `\` on Windows.
    - `/` on its own is the root directory of the whole file system.
    - An absolute path specifies a location from the root of the file system.
    - A relative path specifies a location starting from the current location.
    - `.` on its own means 'the current directory'; `..` means 'the directory above the current one'.

    <Link href="https://swcarpentry.github.io/shell-novice/02-filedir.html">Visit original Carpentries tutorial</Link>

</Alert>

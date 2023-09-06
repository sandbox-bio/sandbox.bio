<script>
import Link from "$components/Link.svelte";
</script>

##### Key Points

- A `for` loop repeats commands once for every thing in a list.
- Every `for` loop needs a variable to refer to the thing it is currently operating on.
- Use `$name` to expand a variable (i.e., get its value). `${name}` can also be used.
- Do not use spaces, quotes, or wildcard characters such as '\*' or '?' in filenames, as it complicates variable expansion.
- Give files consistent names that are easy to match with wildcard patterns to make it easy to select them for looping.
- Use the up-arrow key to scroll up through previous commands to edit and repeat them.
- Use <kbd>Ctrl</kbd>\+<kbd>R</kbd> to search through the previously entered commands.
- Use `history` to display recent commands, and `![number]` to repeat a command by number.

<Link href="https://swcarpentry.github.io/shell-novice/05-loop.html">Visit original Carpentries tutorial</Link>

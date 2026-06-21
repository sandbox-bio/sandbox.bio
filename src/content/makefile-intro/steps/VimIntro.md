<script>
import Execute from "$components/Execute.svelte";
</script>

During this tutorial, we'll use vim to edit our Makefile, so it's worth a quick
refresher. There are many vim tutorials out there, so we'll keep this simple.

1. Basic file editing
2. Save and close
3. Copy, paste, and undo

First, generate a file with some content that we can work with, then print the file to
the console so we can see what it looks like.

<Execute
command={`echo -e 'one, two, three\\nfour, five, six\\nseven, ate, nine' > example.txt`}
/>

<Execute command="cat example.txt" />

Now let's change the typo "ate" to "eight".

<Execute command="vim example.txt" />

VIM has two primary modes: NORMAL and INSERT. NORMAL mode allows you to use shortcuts to
modify and navigate quickly, and INSERT mode allows you to edit the text like a standard
notepad.

1. Edit the file: Hit "i" to enter INSERT mode, make the edit.
2. Save the file: Hit "Esc" to return to NORMAL mode then type `:w` and hit "Enter".
3. Close the editor: Type `:q` and "Enter".

> Tip: You can also save and close with one command, `:wq`, or close without saving
> `:q!`.

Now open the editor back up, and let's paste some text. First, let's try to paste
something from outside of our editor:

1. Copy some text from this tutorial.
2. Hit "i" to enter INSERT mode, then CMD+v (or CTRL+v) to paste.
3. Hit "Esc" to leave INSERT mode, then hit "u" to undo the copy and paste.

This works like you would expect, but you can't copy text within the editor using CMD+c.
But VIM has several shortcuts you can use to copy and paste within the editor.

First, you can copy or cut whole lines:

1. In NORMAL mode, type `yy` to copy the current line.
2. Type `p` to paste it below the cursor (or `P` to paste it above).
3. Delete the current line by typing `dd`.

You can also highlight multiple lines by using the mode VISUAL LINE:

1. Highlight an entire line: `V` (shift+v). (This enables VISUAL LINE mode)
2. Highligh lines above or below by moving the cursor up or down.
3. Hit "Esc" to exit VISUAL LINE.

Or highlight portions of the text by using VISUAL mode:

1. Type `v` to enter VISUAL MODE.
2. Move the cursor around to highlight sections of the text.
3. Exit VISUAL mode by typing either `v` or "Esc".

While text is highlighted, you can copy or cut the selection:

1. Type `v` to enter VISUAL MODE then `e` to highlight to the end of the word.
2. Type `y` to copy or `d` to cut -- this will exit VISUAL mode.
3. Type `p` to paste text after the cursor or `P` (shift+P) to paste before.

Here's a recap of what we covered:

> **Modes**
>
> - `i`: Enter INSERT mode before the cursor (`a` enters after the cursor)
> - `v`: Enter VISUAL mode to highlight
> - `V`: Enter VISUAL LINE mode to highlight lines
> - `Esc`: Leave current mode and return to NORMAL mode
>
> **Commands**
>
> - `:q`: Quit the editor
> - `:w`: Write the buffer to the file (i.e., save the changes)
> - `:q!`: Close without saving
>
> **Shortcuts**
>
> - `u`: Undo the last action (`U` toggles the last action)
> - `yy`: Copy an entire line
> - `dd`: Cut an entire line
> - `y`: In VISUAL mode, copy highlighted text
> - `d`: In VISUAL mode, cut highlighted text
> - `p`: Paste after the cursor (`P` to paste before)

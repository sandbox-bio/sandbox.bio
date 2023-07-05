<script>
import Quiz from "$components/Quiz.svelte";
</script>

We have already mentioned that the Shell lived in a dedicated interface: **the terminal**.

The terminal generally appears as a simple windowed application. Look at the terminal on the right panel.
It is a black frame displaying a **prompt** that usually ends by a `$` sign.
Just after this prompt you will have to type commands.

Several shells can be used to interact within the terminal (e.g. `bash`, `ksh`, `csh`, `zsh`, ...).
The most popular shell is Bash which is the default shell on most Linux systems.
In this course, we will use the Bash shell.

<Quiz id="q1" choices={[
{ valid: true, value: "Shell"},
{ valid: false, value: "Terminal"},
]}>
<span slot="prompt">
Is Bash a shell or a terminal?
</span>
</Quiz>

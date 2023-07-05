<script>
import Quiz from "$components/Quiz.svelte";
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

## The man command

The `man` command is used to get the manual for any command of the Bash shell.

It takes a command name as an argument and opens the manual on the terminal.

Lets try and type:

<Execute command="man ls" />

This manual contains several sections (e.g. NAME, SYNOPSIS, DESCRIPTION).
As indicated in the NAME section this command is used to list a directory content.

```
NAME
       ls - list directory contents
```

The SYNOPSIS section contains the general way of using the command:

```
SYNOPSIS
       ls [OPTION]... [FILE]...
```

The square brackets indicate that both OPTION and FILE are optional for the `ls` command.
Indeed when no options are provided, the `ls` command will simply display the names of files and directories in the current directory without providing the user with additional information (size, owner, creation date...).

And the DESCRIPTION section explains all the possible options of the command.

<Quiz id="q1" choices={[
{ valid: false, value: "list directory contents"},
{ valid: true, value: "use a long listing format"},
{ valid: false, value: "print literal entry names"},
]}>
<span slot="prompt">
What is the meaning of the option `-l` of the `ls` command?
</span>
</Quiz>

## The help option

An other way of getting help is to use the `--help` option after a command name (or, sometimes, `-h` or `-help` or `help`) .
Example, type the following command:

<Execute command="ls --help" />

## Ask Internet

You will also usefull ressources on the Internet:

- <Link href="https://explainshell.com">Explain shell command</Link>
- <Link href="https://stackoverflow.com">A shell forum</Link>

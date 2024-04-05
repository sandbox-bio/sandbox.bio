<script>
import Quiz from "$components/Quiz.svelte";
import Image from "$components/Image.svelte";
</script>

In a Unix operating system, files are arranged in a tree-like structure. In this structure, directories can be seen as branches and files (or empty directories) as leaves. Each file has a unique _path_ in the tree-like structure when starting from the _root_.

Files and directories are accessed through their paths.  
In a path, each successive directory name is separated by a `/`.  
The root of the tree structure is also represented by the first `/` in the path.

There are 2 ways to describe paths: **absolute** and **relative**.

<Image src="/data/ifb-linux-basics-2/absolute_and_relative_paths.png" alt="Absolute and relative paths" />

## Absolute path

An absolute path is described from the root of the tree (ie. beginning with a `/`).
This path is composed of all the names of the different directories from the root of the tree to the target file or directory.

## Relative path and the working directory

With a relative path, you refer to a file or a directory relatively to the directory where you are currently in. We call this directory the _current working directory_.  
A relative path starts from this current working directory, and gives the path from this directory to the target file or directory.

The path of the current working directory can be obtained using the `pwd` command (that stands for **p**rint **w**orking **d**irectory).

Now, type the `pwd` command in the right panel and press <kbd>Enter</kbd>.

<Quiz id="qpwd" choices={[
{ valid: true, value: "An absolute path"},
{ valid: false, value: "A relative path"},
]}>
<span slot="prompt">
What does the `pwd` command return?
</span>
</Quiz>

Look! It begins by a `/`

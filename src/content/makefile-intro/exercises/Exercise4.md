<script>
// Solution:
//    make write-madlib
//    make read-madlib

import Exercise from "$components/Exercise.svelte";

let criteria = [
{
 name: "File <code>word-adverb.txt</code> exists",
 checks: [{
  type: "file",
  path: "word-adverb.txt",
  action: "exists"
 }]
},
{
 name: "File <code>word-noun.txt</code> exists",
 checks: [{
  type: "file",
  path: "word-noun.txt",
  action: "exists"
 }]
},
{
 name: "File <code>word-verb.txt</code> exists",
 checks: [{
  type: "file",
  path: "word-verb.txt",
  action: "exists"
 }]
},
{
 name: "File <code>madlib-beginning.txt</code> contains the story line with the madlb",
 checks: [{
  type: "file",
  path: "madlib-beginning.txt",
  action: "exists",
 }]
},{
 name: "File <code>madlib-middle.txt</code> contains the story line with the madlb",
 checks: [{
  type: "file",
  path: "madlib-middle.txt",
  action: "exists",
 }]
},{
 name: "File <code>madlib-end.txt</code> contains the story line with the madlb",
 checks: [{
  type: "file",
  path: "madlib-end.txt",
  action: "exists",
 }]
}];
</script>

We can upgrade our story to a mad lib! Read through the following data files:

- `Makefile.madlib`
- `generate_word.sh`
- `generate_madlib_part.sh`

Then figure out the targets to build to (a) build the madlib and (b) print out the
story.

<Exercise {criteria} />

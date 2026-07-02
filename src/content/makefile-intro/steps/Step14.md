<!-- TODO: Aside: Make Version and Unused Variables -->

Why was the created file empty in our last step? Did you catch the mistake?

We used upper case for the variable name and lower case for the filename, so `$($*)`
will evaluate to `end` instead of `END`. We'll fix this in the next step, but there's a
couple of small things to consider:

**If you're using `make` v4.4, there's a builtin function that we could use:
`$(upper VALUE)`**

- We're going to do it the hard way as v4.4 was only released a _relatively_ short time
  ago in Oct 2022. Many systems won't have the newer `make` (including the author's
  MacOS). Not a blocker, but we'll use this also introduce a couple of new concepts.

**We don't need `$(STORY_FILES)` anymore.**

- With the original implementation, we used `wildcard` to find the story files on the
  filesystem, and we `patsubst` to build the phony `STORY_RULES` targets. While this
  works, we _could_ end up generating several targets that we don't really need.
- Instead, we explicitly defined the set of targets in `$(STORY_PARTS)` and used the
  pattern substitution shortcut to define both `STORY_FILES` and `STORY_RULES`.
- Because we're using an implicit target `story-%.txt:`, we don't actually need
  `STORY_FILES` any more -- feel free to delete it or use it to build `STORY_RULES`.

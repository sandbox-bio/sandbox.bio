Congratulations! You finished this tutorial!

Let's run through our makefiles to recap what we learned. Your makefiles should look
something like this, but it may not match exactly

###### `Makefile`

```Makefile
# Loads targets from a different makefile
include Makefile.story

# Sets the default target for using `make` without any targets
.DEFAULT_GOAL := story

# Indicates that the targets are _not_ files
.PHONY: hello-world
hello-world:
	@echo "Hello, world."
```

###### `Makefile.story`

```Makefile
# Define the story
BEGINNING := "Once upon a time, there was a chicken."
MIDDLE := "It crossed the road."
END := "And lived happily ever after."

# Build the story targets
STORY_PARTS := beginning middle end
STORY_FILES := $(STORY_PARTS:%=story-%.txt)
STORY_RULES := $(patsubst %.txt,%,$(STORY_FILES))

to_upper = $(shell echo "$(1)" | tr '[:lower:]' '[:upper:]')

.PHONY: story $(STORY_RULES)

# debug rule
check:
	@echo "parts: ${STORY_PARTS}"
	@echo "files: ${STORY_FILES}"
	@echo "rules: ${STORY_RULES}"

clean: clean-story
clean-story:
	-@rm -f story-*.txt

# Build the story files
story-%.txt:
	@echo "$($(call to_upper,$*))" > $@

# Phony rules to print the story to the CLI
$(STORY_RULES): story-%: story-%.txt
	@cat "story-${*}.txt"

# Ensure the story is written out in order
story-middle: story-beginning
story-end: story-middle

# Print the story to the command line
story: story-end
	@echo "The End."
```

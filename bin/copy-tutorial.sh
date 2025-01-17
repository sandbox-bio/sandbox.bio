#!/bin/bash

# Script to copy a tutorial from another repo
# Make sure to run from the "app/" folder.
# Usage: ./bin/copy-tutorial.sh /source/path/ destination-tutorial-name

# Example:
# ./bin/copy-tutorial.sh "/Users/robert/Documents/dev/sandboxbioscenarios/linux_basics_session01" "ifb-linux-basics-1" "false"
# ./bin/copy-tutorial.sh "/Users/robert/Documents/dev/sandboxbioscenarios/linux_basics_session02" "ifb-linux-basics-2" "false"
# ./bin/copy-tutorial.sh "/Users/robert/Documents/dev/sandboxbioscenarios/linux_basics_session03" "ifb-linux-basics-3" "false"
# ./bin/copy-tutorial.sh "/Users/robert/Documents/dev/sandboxbioscenarios/linux_basics_session04" "ifb-linux-basics-4" "false"

# Update an existing tutorial from a given file path (e.g. clone repo in that path)
TUTORIAL_SOURCE_DIR=${1?Missing source tutorial path}
TUTORIAL_DEST_ID=${2?Missing destination tutorial name}
TUTORIAL_LISTED=${3:-"true"}

TUTORIAL_SOURCE_ID=$(basename ${TUTORIAL_SOURCE_DIR})
TUTORIAL_DEST_DIR="src/content/${TUTORIAL_DEST_ID}"

# Sync file contents
mkdir -p $TUTORIAL_DEST_DIR
cp -r ${TUTORIAL_SOURCE_DIR}/* src/content/$TUTORIAL_DEST_ID

# Update config.js to use tutorial ID of interest
sed -i.bak "s/${TUTORIAL_SOURCE_ID}/${TUTORIAL_DEST_ID}/g" $TUTORIAL_DEST_DIR/config.js
sed -i.bak "s/listed: .*,/listed: ${TUTORIAL_LISTED},/" $TUTORIAL_DEST_DIR/config.js

# Don't need prefix "data/tutorial-name" under `files`
sed -i.bak "s|data/${TUTORIAL_DEST_ID}/||" $TUTORIAL_DEST_DIR/config.js

# Fix import Quiz from "components/Quiz.svelte"; ---> $components/
sed -i.bak 's|from "components/|from "$components/|' $TUTORIAL_DEST_DIR/*/*.md

# Fix yellow belt/orange belt
sed -i.bak 's/difficulty:.*/difficulty: ["beginner"],/' $TUTORIAL_DEST_DIR/config.js

# Update tutorials.js
importName=${TUTORIAL_DEST_ID//-}
importStatement="/${TUTORIAL_DEST_ID}/config.js\";"  # partial match because v1 tutorials used /tutorials/, not $content/
# Partial match because v1 tutorials
if [[ $(grep -c "$importStatement" src/stores/tutorials.js) -lt 1 ]]; then
    # Include import statement
    sed -i.bak '2i\'$'\n'''"import { config as ${importName} } from \"\$content$importStatement"''$'\n' src/stores/tutorials.js

    # Add to list of tutorials
    sed -i.bak "s|\t// All tutorials|\t$importName,\n\t// All tutorials|" src/stores/tutorials.js
    rm src/stores/tutorials.js.bak
fi

./bin/build.sh

rm $TUTORIAL_DEST_DIR/config.js.bak
rm src/stores/tutorials.js.bak
rm $TUTORIAL_DEST_DIR/*/*.md.bak

npm run format

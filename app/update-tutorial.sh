#!/bin/bash

# Script to update a tutorial
# Make sure to run from the "app/" folder.
# Usage: ./update-tutorial.sh /source/path/ destination-tutorial-name

# Update an existing tutorial from a given file path (e.g. clone repo in that path)
TUTORIAL_SOURCE_DIR=${1?Missing source tutorial path}
TUTORIAL_DEST_ID=${2?Missing destination tutorial name}
TUTORIAL_LISTED=${3:-"true"}

TUTORIAL_SOURCE_ID=$(basename ${TUTORIAL_SOURCE_DIR})
TUTORIAL_DEST_DIR="src/tutorials/${TUTORIAL_DEST_ID}"

# Sync file contents
mkdir -p $TUTORIAL_DEST_DIR
cp -r ${TUTORIAL_SOURCE_DIR}/* src/tutorials/$TUTORIAL_DEST_ID

# Update config.js to use tutorial ID of interest
sed -i.bak "s/${TUTORIAL_SOURCE_ID}/${TUTORIAL_DEST_ID}/g" $TUTORIAL_DEST_DIR/config.js
sed -i.bak "s/listed: .*,/listed: ${TUTORIAL_LISTED},/" $TUTORIAL_DEST_DIR/config.js

# Update tutorials.js
importName=${TUTORIAL_DEST_ID//-}
importStatement="import { config as ${importName} } from \"tutorials/${TUTORIAL_DEST_ID}/config.js\";"
if [[ $(grep -c "$importStatement" src/stores/tutorials.js) -lt 1 ]]; then
    # Include import statement
    sed -i.bak '2i\'$'\n'''"$importStatement"''$'\n' src/stores/tutorials.js

    # Add to list of tutorials
    sed -i.bak "s|\t// Add tutorials here|\t$importName,\n\t// Add tutorials here|" src/stores/tutorials.js
fi
./build.sh

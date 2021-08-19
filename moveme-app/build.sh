#!/bin/bash

# Add build date and commit ID
echo "<!-- Build: $(TZ=America/Los_Angeles date); $(git rev-parse HEAD) -->" >> public/index.html

# Copy over data we need to be available via URL for mounting data as URLs, and for IGV.js
mkdir -p public/data
for tutorial in $(ls -d src/tutorials/*);
do
	dest=public/data/$(basename $tutorial)
	mkdir -p $dest
	cp $tutorial/data/* $dest
done

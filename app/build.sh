#!/bin/bash

# Add build date and commit ID
echo "<!-- Build: $(TZ=America/Los_Angeles date); $(git rev-parse HEAD) -->" >> public/index.html

# Create mock routes (only used for local testing; there's a router server-side to handle this)
PAGES=("tutorials" "playground" "rosalind")
for page in ${PAGES[@]}; do
	mkdir -p "public/$page"
	cp public/index.html public/$page/index.html
done

# Copy over data we need to be available via URL for mounting data as URLs, and for IGV.js
mkdir -p public/data
for tutorial in $(ls -d src/tutorials/*/);
do
	[[ ! -e $tutorial/data ]] && continue;
	dest=public/data/$(basename $tutorial)
	mkdir -p $dest
	cp $tutorial/data/* $dest
done

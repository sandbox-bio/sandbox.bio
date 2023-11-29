#!/bin/bash

# Copy over data we need to be available via URL for mounting data as URLs, and for IGV.js
mkdir -p static/data
for tutorial in $(ls -d src/content/*/);
do
	[[ ! -e $tutorial/data ]] && continue;
	dest=static/data/$(basename $tutorial)
	mkdir -p $dest
	cp -r $tutorial/data/* $dest
done

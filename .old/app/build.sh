#!/bin/bash

# Copy over data we need to be available via URL for mounting data as URLs, and for IGV.js
mkdir -p public/data
for tutorial in $(ls -d src/tutorials/*/);
do
	[[ ! -e $tutorial/data ]] && continue;
	dest=public/data/$(basename $tutorial)
	mkdir -p $dest
	cp -r $tutorial/data/* $dest
done

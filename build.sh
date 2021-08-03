#!/bin/bash

# Crude way to support routing: copy index.html and 
mkdir -p public/tutorials public/playground
cp public/index.html public/tutorials/
cp public/index.html public/playground/

# Copy over data we need to be available via URL for mounting data as URLs, and for IGV.js
mkdir -p public/data

for tutorial in $(ls -d src/tutorials/*);
do
	dest=public/data/$(basename $tutorial)
	mkdir -p $dest
	cp $tutorial/data/* $dest
done

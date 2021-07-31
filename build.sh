#!/bin/bash

# Crude way to support routing: copy index.html and 
mkdir -p public/tutorials public/playground
cp public/index.html public/tutorials/
cp public/index.html public/playground/

# Copy over data we need to be available via URL for IGV.js
mkdir -p public/data
cp src/tutorials/bedtools-intro/data/{cpg,exons,gwas,hesc.chromHmm}.bed public/data/

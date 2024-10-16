#!/bin/bash

# This script simulates the behavior of curl for a specific URL pattern.
# It accepts an optional -s argument to suppress stderr output.

# Check if -s option is present
if [[ $1 == "-s" ]]; then
    # If -s option is present, shift all arguments to the left
    shift
    # Set stderr to /dev/null to suppress output
    exec 2>/dev/null
fi

echo "Simulating downloading file..." >&2
sleep 1

# Store the URL argument in the variable url
url=$1

# Check if the URL matches the pattern for ChEBI IDs
if [[ $url == https://www.ebi.ac.uk/chebi/* ]]; then
    # Extract the ChEBI ID from the URL using a regular expression
    chebiId=$(echo $url | grep -oP '(?<=chebiId=)[^&]*')
    # Print the contents of the corresponding CSV file
    /usr/local/bin/curl.ori "/data/datepro-02-data-retrieval/chebi_${chebiId}_xrefs_UniProt.csv"
else
    # If the URL does not match the expected pattern, print an error message
    echo "URL for an unknown location"
fi

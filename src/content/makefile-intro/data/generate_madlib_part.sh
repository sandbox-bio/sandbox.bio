#!/bin/bash

MADLIB_PART=$1
MADLIB_FILE=madlib-${MADLIB_PART}.txt

if [ $# -lt 1 ]; then
    echo "Error: Missing required argument"
    echo "Usage: $0 MADLIB_PART"
    exit 1
fi

if [[ "${MADLIB_PART}" == "beginning" ]]; then
    word=$(cat word-noun.txt)
    echo "Once upon a time, there was a ${word}." > ${MADLIB_FILE}
elif [[ "${MADLIB_PART}" == "middle" ]]; then
    word=$(cat word-veb.txt)
    echo "It ${word}ed the road." > ${MADLIB_FILE}
elif [[ "${MADLIB_PART}" == "end" ]]; then
    word=$(cat word-adverb.txt)
    echo "And lived ${word} on the other side." > ${MADLIB_FILE}
else
    echo "Unknown story part: ${MADLIB_PART}"
    exit 1
fi
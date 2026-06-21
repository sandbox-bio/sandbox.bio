#!/bin/bash

PART_OF_SPEECH=$1
POS_FILE=word-${PART_OF_SPEECH}.txt

if [ $# -lt 1 ]; then
    echo "Error: Missing required argument"
    echo "Usage: $0 PART_OF_SPEECH"
    exit 1
fi

if [ -f ${POS_FILE} ]; then
    # Allow use of current value
    CURRENT=$(cat ${POS_FILE})
    read -r -p "Enter a ${PART_OF_SPEECH} (or nothing to reuse current value ${CURRENT}): " user_value
    if [[ -z ${user_value} ]]; then
        # Do nothing -- leave as is
        echo "Reusing current value: ${PART_OF_SPEECH}=${CURRENT}"
        exit 0
    fi
else
    # Get new value
    read -r -p "Enter a ${PART_OF_SPEECH}: " user_value
fi

echo "Using given value: ${PART_OF_SPEECH}=${user_value}"
echo "$user_value" > "${POS_FILE}"
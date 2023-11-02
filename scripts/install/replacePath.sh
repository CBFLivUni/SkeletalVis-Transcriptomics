#!/bin/bash

# Check for the correct number of command-line arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input_file> <search_path> <replacement_path>"
  exit 1
fi

input_file="$1"
search_path="$2"
replacement_path="$3"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Error: The input file '$input_file' does not exist."
  exit 1
fi

# Perform the search and replace using sed
sed -i "s|$search_path|$replacement_path|g" "$input_file"

echo "Replacement complete. $search_path has been replaced with $replacement_path in $input_file."

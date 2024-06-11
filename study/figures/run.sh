#!/bin/bash

# Find all .py files recursively from the current directory
find "$(pwd)" -name "*.py" | while read -r file; do
  echo "Working on $file"

  # Get the absolute directory containing the file
  dir=$(dirname "$file")
  
  # Move to the directory
  cd "$dir" || exit
  
  # Execute the Python script
  python "$(basename "$file")"
  
  # Move back to the original directory
  cd - > /dev/null || exit
done

#!/bin/bash

# Check if 'renv.lock' file exists and remove it if it does
if [ -f renv.lock ]; then
    echo "'renv.lock' file found. Removing..."
    rm -f -- renv.lock
else
    echo "'renv.lock' file not found. Skipping..."
fi

# Define the target folder
target_folder="src"

# Find and delete files with specific extensions in the target folder
echo "Searching for files with extensions '.o', '.dll', or '.so' in $target_folder..."
find "$target_folder" -type f \( -name "*.o" -o -name "*.dll" -o -name "*.so" \) -exec rm -f {} \;
echo "Deletion complete."

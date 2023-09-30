#!/bin/bash

target_folder="../src"

find "$target_folder" -type f \( -name "*.o" -o -name "*.dll" -o -name "*.so" \) -exec rm -f {} \;

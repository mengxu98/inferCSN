#!/bin/bash

# Check the status of your local repository
git status

# Add the changes to the staging area
git add .

# Commit the changes with a descriptive message
git commit -m "commit message"

# Switch to the dev branch
git checkout dev

# Pull the latest changes from the remote dev branch
git pull origin dev

# Push your local changes to the remote dev branch
git push origin dev

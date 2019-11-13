#!/bin/bash

#
# Running this script will build and push all documentation
# to GitHub Pages.
# 

GIT_BRANCH=$(git branch | grep \* | cut -d ' ' -f2)

git symbolic-ref HEAD refs/heads/gh-pages
cd ..
mkdocs build
rm .git/index
git add html/*
git mv html/* .
touch .nojekyll && git add .nojekyll
git commit -m "Update documentation"
git push origin gh-pages
git checkout -f ${GIT_BRANCH}
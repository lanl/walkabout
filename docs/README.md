# Documentation - Walkabout (V1)
### Approved for public release; LA-UR-11-01952

Repository is ready for development. This version starts from a debugged version of Soumi’s code containing OpenMP parallelization (See changelog for further detail.) <br />

Begin date: 07/12/16 <br />
End date: 
Locations: https://github.com/lanl/walkabout.git; /scratch/ymp/smckinney/walkabout/walkabout/ <br />

__Update in Progress__
[Scott Painter's User's Manual](./WalkaboutUM.pdf)

### Developer's Guide to Compiling Documentation

Install MkDocs with:

```
pip install mkdocs
```

Install `mdx_math` extension with:

```
pip install python-markdown-math
```

Install `material` theme with:

```
pip install mkdocs-material
```

From the top level of this directory (`..`), run:

```
mkdocs build
```

to compile HTML documentation, or

```
mkdocs serve
```

to run a local server (useful for development, as changes are updated
automatically).

## Pushing to GitHub

To appear on GitHub pages, HTML must be built on the `gh-pages` branch.

Run the following command exactly from the `walkabout/` root directory to
update and push documentation:

```
git checkout gh-pages && git merge master \
  && rm -Rf about assets user-guide images search supplementary-docs \
  && mkdocs build && mv html/* . \
  && git rm -r test/test4 \
  && git add . && git commit -m "Update html" \
  && git push origin gh-pages \
  && git checkout master
```

*Note: remember to commit your changes on `master` before running that command!*
# Building and reviewing the documentation

Use a supported Python environment with pele and its native dependencies
installed (see the repository README). The docs build imports the current
package, including compiled extensions; it does not mock pele. Install the docs
and GUI extras and Doxygen. On Ubuntu the OpenGL imports also need GLU and GLUT:

```sh
sudo apt-get install doxygen rsync libglu1-mesa freeglut3-dev
python -m pip install '.[docs,gui]' \
  -Csetup-args=-Dgui=enabled -Csetup-args=-Dlammps=disabled -Ccompile-args=-j1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen \
  make -C doc check SPHINXOPTS='-n -W --keep-going -E -a'
python -I -m unittest discover -s doc/tests
```



`doc/Makefile` runs Doxygen first using `doc/doxygen/Doxyfile`. It extracts the
C++ sources into `doc/doxygen/_build/xml` and a full HTML reference. Sphinx's
Breathe extension reads that XML and renders the core interfaces in
`cpp_reference.rst`, with Sphinx search, navigation, and C++ reference targets.
The complete Doxygen HTML is copied under `doc/_build/html/doxygen`; its
navigation links back to Sphinx. References to types outside the selected
Breathe classes resolve to that complete reference using Doxygen's XML index.

For Python, Sphinx autodoc imports the modules and reads their docstrings.
Autosummary follows the explicit catalogs in the `.rst` pages and package
docstrings, using `doc/_templates/autosummary` to generate class, method, and
attribute pages under `doc/generated`. New exports need an autosummary entry
to receive a reference page. Numpydoc renders the NumPy-style docstrings.
The version comes from `pyproject.toml`. Generated sources, Doxygen output,
and final HTML are ignored by Git; edit their inputs instead.

`make -C doc check` builds the combined site, verifies representative Python
pages and actual Breathe C++ targets, then checks local HTML links, anchors,
images, scripts, and stylesheets. `-W --keep-going` makes Sphinx warnings fail
the build; `-n` also checks unresolved API references. Doxygen's existing source-comment warnings are reported separately;
they are not fatal. External URLs can be checked with `make -C doc linkcheck`,
but that requires network access and can fail for unavailable remote sites.

For a clean regeneration, use `make -C doc clean` before the build. To review
the result, serve `doc/_build/html` locally or download the `pele-documentation`
artifact from the Documentation workflow. Artifacts are retained for 14 days.

## CI and publication

`.github/workflows/docs.yml` builds and checks documentation for pushes to
`docs-autogeneration` and `main`, PRs targeting `main`, and manual runs. Every
successful build uploads both the review site and a deployable Pages artifact.
PRs cannot deploy. CI builds the current checkout with native dependencies and
serial compilation.

Publication is disabled unless the repository variable `DOCS_PUBLISH_ENABLED`
is exactly `true`. Pages must use **GitHub Actions** as its publishing source.
Successful `main` builds then refresh the website; manual main runs use the
same deployment path. No custom deployment secret is needed: the deploy job
uses GitHub's Pages and OIDC permissions.

For an owner-approved first deployment before merging, temporarily set the
repository variable `DOCS_PUBLISH_BRANCH=docs-autogeneration` as well, then push
to that branch. Remove that variable after the deployment succeeds so later
docs-branch edits only build and check. Keep the old `gh-pages` branch for
rollback. Artifacts expire after 14 days; rerun the build if necessary.

Publishing is a separate owner-approved action from preparing a PR. Repository
and environment review requirements can be configured separately; the
`github-pages` environment name in YAML alone does not require approval.
Installation and getting-started pages are the first website review priority;
the older custom-potential and system-class tutorials still need a separate
Python 3 content review.

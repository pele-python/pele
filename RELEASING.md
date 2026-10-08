# Releasing pele

Release tested commits from `main` using matching version tags and GitHub Releases.
Publish the source distribution and repaired wheels to PyPI. Wheels support
CPython 3.11–3.14 on Linux x86_64 (manylinux_2_28) and macOS 14+ on Intel and
Apple Silicon. Windows is not supported. Conda-forge packages follow recipe review.

## Prepare a version

Keep the version identical in `pyproject.toml`, `meson.build`, and
`conda-recipe/recipe.yaml`. The initial version is `0.1.1`.

Run both the Tests and Release artifacts workflows successfully on the release
commit. Release artifacts builds the source archive once, then builds and tests
all 12 wheels from that archive on Linux and macOS with Python 3.11–3.14,
including generated GUI forms. Wheel jobs build pinned native dependencies,
including macOS OpenMP. macOS jobs reuse completed dependency builds only when
the source recipe, architecture, compilers, SDK, OS and deployment target match.
macOS jobs install SHA-verified Sonoma GCC bottles
to keep the bundled Fortran runtimes compatible with macOS 14. Jobs repair the wheels
with auditwheel/delocate, then run installed-package tests with native-library
search paths cleared. The tests reject libraries loaded from the build prefix or
Homebrew. All three wheel jobs must pass before publication. Distribution builds
use `native=false`; use `-Dnative=true` only for local CPU-specific builds.

The native wheel dependency builder is `ci/build-wheel-deps.sh`; its source URLs
and SHA256 values are pinned. `LICENSE-THIRD-PARTY` records redistributed-library
notices. Before changing dependencies, verify their licenses, source availability,
and ABI/macOS deployment requirements. Delocate must retain its deployment-target
checks; auditwheel must successfully repair to the declared manylinux policy.

To validate locally, first install C/C++/Fortran compilers, SUNDIALS, Eigen, LAPACK
headers, and OpenMP, then run:

```sh
python -m pip install build twine
python -m build -Csetup-args=-Dnative=false -Csetup-args=-Dlammps=disabled
python -m twine check --strict dist/pele-0.1.1.tar.gz
```

The default `python -m build` builds a source archive and then builds its wheel
from that archive. Keep the checkout committed: Meson creates the source archive
from Git's committed files. The compiled LAMMPS integration is excluded from the
generic distribution; the Python LAMMPS fallback remains available through the
`lammps` extra.

For a local Linux conda build outside conda-forge, supply its standard-library
variants explicitly (conda-forge injects these during staging):

```sh
rattler-build build -r conda-recipe/recipe.yaml -c conda-forge \
  --variant c_stdlib=sysroot --variant c_stdlib_version=2.17
```

## Configure publishing accounts

Create accounts on [PyPI](https://pypi.org/) and [TestPyPI](https://test.pypi.org/).
Configure a pending Trusted Publisher on each service with these values:

| Field | PyPI | TestPyPI |
| --- | --- | --- |
| Project name | `pele` | `pele` |
| GitHub owner | `pele-python` | `pele-python` |
| Repository | `pele` | `pele` |
| Workflow filename | `release.yml` | `release.yml` |
| Environment | `pypi` | `testpypi` |

Create the corresponding GitHub environments under Settings → Environments if
they do not exist. Trusted Publishing uses GitHub's identity token; no PyPI API
key is needed. A pending publisher does not reserve a package name until its
first upload.

## Test the upload

On GitHub, open Actions → Release artifacts → Run workflow. Select `main` and
check the TestPyPI input. Publishing begins only after every validation job passes.
A dispatch with the input unchecked validates artifacts without publishing.

Test a binary installation in a fresh supported environment:

```sh
python -m pip install --only-binary=:all: --no-deps --index-url https://test.pypi.org/simple/ \
  pele==0.1.1
```

Install the runtime dependencies listed in `pyproject.toml` before using
`--no-deps`. Use a new version for each revised TestPyPI upload; distribution
filenames cannot be replaced.

## Publish the release

Fetch `main` from `pele-python/pele` into a clean release checkout. Once the
workflows pass and publishing accounts are configured, create and push an
annotated tag:

```sh
git tag -a v0.1.1 -m 'pele 0.1.1'
git push https://github.com/pele-python/pele.git v0.1.1
```

Create a GitHub Release for `v0.1.1` and publish it. The Release artifacts workflow
requires the tag to match the package version, validates the distribution, and
publishes the source archive and all repaired wheels to PyPI. Creating a tag
alone does not upload a package.

## Submit the conda-forge recipe

The recipe in this repository uses a local source path for development builds.
After the PyPI source archive is published, find its exact URL and SHA256 at
`https://pypi.org/pypi/pele/0.1.1/json` (the `urls` entry whose `packagetype` is
`sdist`). Replace the recipe's entire `source` block with:

```yaml
source:
  url: <the published sdist URL>
  sha256: <its SHA256 digest>
```

Copy the resulting recipe to `recipes/pele/recipe.yaml` in a fork of
[conda-forge/staged-recipes](https://github.com/conda-forge/staged-recipes), and
submit a pull request. Keep `spraharsh` under `recipe-maintainers`. The source URL
and checksum refer to the same tested distribution used for PyPI. Conda-forge
review and feedstock builds are required before `conda install -c conda-forge
pele` becomes available.

For subsequent releases, publish a new version and update the version, source
checksum, and dependencies in `conda-forge/pele-feedstock`. Start a maintenance
branch only when an older supported release needs fixes alongside newer work.

[PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/)
and [conda-forge submission instructions](https://conda-forge.org/docs/maintainer/adding_pkgs/)
describe the account and review steps.

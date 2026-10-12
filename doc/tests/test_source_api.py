import importlib.util
import json
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from docutils import nodes
from docutils.core import publish_doctree


class SourceApiTests(unittest.TestCase):
    def generator(self):
        path = Path(__file__).resolve().parents[1] / "sphinxext/source_api.py"
        self.assertTrue(path.is_file(), "The source API generator has not been implemented")
        spec = importlib.util.spec_from_file_location("source_api", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_public_python_apis_in_private_modules_exclude_private_and_local_helpers(self):
        generator = self.generator()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "pele/_example.py"
            source.parent.mkdir()
            source.write_text('''"""Library overview.

.. automodule:: unavailable_backend
   :members:
"""
import unavailable_backend
raise RuntimeError("Documentation must not import this module")

class Example:
    """A documented class."""
    def __init__(self):
        """A documented public constructor."""
    def __call__(self):
        """A documented public call protocol."""
    def __getitem__(self, key):
        """A documented public indexing protocol."""
    def method(self, x=2):
        """A documented public method."""
    def _hidden(self, x=2):
        """A documented private method."""
    def __hidden(self):
        """A documented name-mangled private method."""

class _Private:
    """A private class must not be published."""
    def public_name(self):
        """A method of a private class must not be published."""

def outer(value):
    """A documented function."""
    class Nested:
        """A documented nested class."""
    def local():
        """A documented nested function."""
''')
            skipped = root / "pele/tests/test_example.py"
            skipped.parent.mkdir()
            skipped.write_text('def test_example():\n    """Do not publish tests."""\n')
            (root / "pele/untracked.py").write_text('def hidden():\n    """Do not publish untracked code."""\n')
            private_only = root / "pele/_private_only.py"
            private_only.write_text('"""Private-only module prose."""\ndef _helper():\n    """A private helper."""\n')
            public_prose = root / "pele/intro.py"
            public_prose.write_text('"""Public module prose."""\n')
            output = root / "doc/generated/source_api"
            with patch.object(generator, "tracked_sources", return_value=[source, skipped, private_only, public_prose]):
                manifest = generator.generate(root, output, revision="test-revision")
                second = generator.generate(root, output, revision="test-revision")
            self.assertEqual(manifest, second)
            self.assertEqual({entry["fullname"] for entry in manifest}, {
                "pele._example", "pele._example.Example", "pele._example.Example.__init__",
                "pele._example.Example.__call__", "pele._example.Example.method",
                "pele._example.Example.__getitem__",
                "pele._example.outer", "pele.intro",
            })
            page = (output / "pele._example.rst").read_text()
            self.assertNotIn(".. automodule::", page)
            self.assertIn("Library overview.", page)
            self.assertIn("A documented public method.", page)
            self.assertNotIn("A documented private method.", page)
            self.assertNotIn("name-mangled private method", page)
            self.assertNotIn("A documented nested", page)
            self.assertFalse((output / "pele._private_only.rst").exists())
            self.assertIn("blob/test-revision/pele/_example.py#L", page)
            links = "\n".join(line.strip() for line in page.splitlines() if line.lstrip().startswith("`Source <"))
            tree = publish_doctree(links, settings_overrides={"warning_stream": StringIO()})
            self.assertFalse([message.astext() for message in tree.findall(nodes.system_message) if message["level"] >= 2])
            for entry in manifest:
                symbol_page = output / (entry["docname"].rsplit("/", 1)[-1] + ".rst")
                self.assertIn(f".. _{entry['anchor']}:", symbol_page.read_text())
            self.assertEqual(json.loads((output / "manifest.json").read_text()), manifest)

    def test_public_cython_types_in_private_modules_exclude_private_scopes(self):
        generator = self.generator()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "pele/_native.pyx"
            source.parent.mkdir()
            source.write_text('''import unavailable_backend
cimport unavailable_backend
raise RuntimeError("Do not execute Cython module code")

cdef class Native:
    """A documented extension class."""
    @unknown_decorator
    def __cinit__(self, double value=1):
        """A documented constructor."""
    cdef double _hidden(self, double value) noexcept:
        """A documented native method."""
    cdef double private_method(self, double value) noexcept:
        """An implicitly private native method."""
    cpdef double public_method(self, double value) noexcept:
        """A documented public native method."""

cdef class _HiddenNative:
    """A private extension class must not be published."""
    def public_name(self):
        """A method of a private extension class must not be published."""

cpdef double bridge(
    double value,
    int count=2,
):
    """A documented bridge."""
    class Nested:
        """A nested Cython class must not be published."""

cdef public double public_export(double value):
    """A documented exported native function."""

cdef api double api_export(double value):
    """A documented Cython API function."""
''')
            declaration = root / "pele/helpers.pxd"
            declaration.write_text('''cdef inline double helper(double x) noexcept:
    """A documented declaration helper."""
    return x
''')
            output = root / "doc/generated/source_api"
            with patch.object(generator, "tracked_sources", return_value=[source, declaration]):
                manifest = generator.generate(root, output, revision="test-revision")
            self.assertEqual({entry["fullname"] for entry in manifest}, {
                "pele._native.Native", "pele._native.Native.__cinit__",
                "pele._native.Native.public_method", "pele._native.bridge",
                "pele._native.public_export", "pele._native.api_export",
            })
            page = (output / "pele._native.rst").read_text()
            self.assertIn("A documented constructor.", page)
            self.assertIn("double value", page)
            self.assertIn("noexcept", page)
            self.assertNotIn("A documented native method.", page)
            self.assertNotIn("An implicitly private native method.", page)
            self.assertNotIn("A nested Cython class", page)
            self.assertFalse((output / "pele.helpers.rst").exists())

    def test_module_catalog_removal_preserves_heading_levels_and_spacing(self):
        generator = self.generator()
        doc = """Utilities
+++++++++
Helper functions
================

Rotations
---------
.. automodule:: unavailable_backend
   :members:

Benchmarking
------------
.. autosummary::
   :toctree: generated/

   unavailable_backend.Function
"""
        prose = generator.module_prose(doc)
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "pele/utils.py"
            source.parent.mkdir()
            source.write_text(repr(doc) + "\n")
            output = root / "doc/generated/source_api"
            with patch.object(generator, "tracked_sources", return_value=[source]):
                generator.generate(root, output)
            page = (output / "pele.utils.rst").read_text()
            page = page.replace(".. py:currentmodule:: pele.utils", "")
            tree = publish_doctree(page, settings_overrides={"warning_stream": StringIO()})
            self.assertFalse(list(tree.findall(nodes.system_message)))
        self.assertIn("\n\nBenchmarking", prose)
        self.assertNotIn("unavailable_backend", prose)

    def test_manifest_is_not_copied_when_strict_build_has_warnings(self):
        generator = self.generator()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "generated/source_api"
            source.mkdir(parents=True)
            (source / "manifest.json").write_text("[]\n")
            output = root / "html"
            output.mkdir()
            callbacks = {}
            app = SimpleNamespace(
                srcdir=root, outdir=output, statuscode=0, _warncount=1,
                builder=SimpleNamespace(format="html"),
                add_config_value=lambda *args: None,
                connect=lambda event, callback, **kwargs: callbacks.update({event: callback}),
            )
            generator.setup(app)
            callbacks["build-finished"](app, None)
            self.assertFalse((output / "source_api_manifest.json").exists())
            app._warncount = 0
            callbacks["build-finished"](app, None)
            self.assertEqual((output / "source_api_manifest.json").read_text(), "[]\n")

    def test_unrecognized_docstring_sections_keep_their_content(self):
        generator = self.generator()
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "pele/unusual.py"
            source.parent.mkdir()
            source.write_text('''def unusual(value):
    """A documented function.

    Input
    -----
    A description that Numpydoc does not recognize as a standard section.
    """
''')
            output = root / "doc/generated/source_api"
            with patch.object(generator, "tracked_sources", return_value=[source]):
                generator.generate(root, output)
            self.assertIn("A description that Numpydoc does not recognize", (output / "pele.unusual.rst").read_text())

if __name__ == "__main__":
    unittest.main()

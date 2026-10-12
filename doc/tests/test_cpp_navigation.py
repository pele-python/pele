import importlib.util
from pathlib import Path
import posixpath
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from jinja2 import Template


class CppNavigationTests(unittest.TestCase):
    def extension(self):
        path = Path(__file__).resolve().parents[1] / "sphinxext/cpp_navigation.py"
        spec = importlib.util.spec_from_file_location("cpp_navigation", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def fixture(self, directory):
        source = Path(directory)
        output = source / "generated/cpp"
        output.mkdir(parents=True)
        nodes = [
            ("GradientOptimizer", [], "gradient.rst"),
            ("LBFGS", [("public", "classpele_1_1GradientOptimizer", "pele::GradientOptimizer")], "lbfgs.rst"),
            ("NestedOptimizer", [("protected", "classpele_1_1LBFGS", "pele::LBFGS")], output / "nested/optimizer.rst"),
            ("BasePotential", [], "base_potential.rst"),
            ("Potential", [("private", "classpele_1_1BasePotential", "pele::BasePotential")], "potential.rst"),
            ("CycleA", [("public", "classpele_1_1Potential", "pele::Potential"),
                        ("public", "classpele_1_1CycleB", "pele::CycleB")], "cycle_a.rst"),
            ("CycleB", [("public", "classpele_1_1CycleA", "pele::CycleA")], "cycle_b.rst"),
            ("UnrelatedA", [("public", "classpele_1_1UnrelatedB", "pele::UnrelatedB")], "unrelated_a.rst"),
            ("UnrelatedB", [("public", "classpele_1_1UnrelatedA", "pele::UnrelatedA")], "unrelated_b.rst"),
        ]
        nodes = [SimpleNamespace(kind="class", name="pele::" + name, refid="classpele_1_1" + name,
                                 base_compounds=bases, file_name=str(filename))
                 for name, bases, filename in nodes]
        for node in nodes:
            page = output / node.file_name
            page.parent.mkdir(parents=True, exist_ok=True)
            page.write_text(f"Canonical definition of {node.name}.\n")
        root_page = output / "index.rst"
        root_page.write_text("C++ API reference\n=================\n")
        graph = SimpleNamespace(root_directory=str(output), full_root_file_path=str(root_page),
                                all_nodes=nodes, class_like=[node for node in nodes if node.name != "pele::NestedOptimizer"])
        return SimpleNamespace(srcdir=str(source), exhale_root=graph,
                               config=SimpleNamespace(html_context={"existing": "preserved"}))

    def test_topics_follow_all_public_class_inheritance_and_link_canonical_pages(self):
        extension = self.extension()
        with TemporaryDirectory() as directory:
            app = self.fixture(directory)
            output = Path(app.exhale_root.root_directory)
            definitions = {page: page.read_text() for page in output.rglob("*.rst") if page.name != "index.rst"}
            extension.generate(app)
            extension.generate(app)
            optimizers = (output / "optimizers.rst").read_text()
            potentials = (output / "potentials.rst").read_text()
            self.assertIn(":doc:`pele::GradientOptimizer </generated/cpp/gradient>`", optimizers)
            self.assertIn(":doc:`pele::LBFGS </generated/cpp/lbfgs>`", optimizers)
            self.assertIn(":doc:`pele::NestedOptimizer </generated/cpp/nested/optimizer>`", optimizers)
            self.assertEqual(optimizers.count(":doc:"), 3)
            for name in ("BasePotential", "Potential", "CycleA", "CycleB"):
                self.assertIn("pele::" + name, potentials)
            self.assertEqual(potentials.count(":doc:"), 4)
            self.assertNotIn("Unrelated", potentials)
            self.assertNotIn(".. doxygen", optimizers + potentials)
            self.assertEqual((output / "index.rst").read_text().count("Browse by topic"), 1)
            self.assertEqual(definitions, {page: page.read_text() for page in definitions})
            self.assertEqual(app.config.html_context["existing"], "preserved")

    def test_sidebar_categories_link_correctly_from_every_page_depth(self):
        extension = self.extension()
        with TemporaryDirectory() as directory:
            app = self.fixture(directory)
            extension.generate(app)
            navigation = app.config.html_context["cpp_navigation"]
            self.assertEqual([(entry["label"], entry["anchor"]) for entry in navigation], [
                ("Overview", ""), ("Optimizers", ""), ("Potentials", ""),
                ("Classes & interfaces", "class-hierarchy"), ("Namespaces", "namespaces"),
                ("Functions", "functions"), ("Source files", "file-hierarchy"),
            ])
            template = Template((Path(__file__).resolve().parents[1] / "_templates/api_navigation.html").read_text())
            for pagename in ("index", "generated/cpp/index", "generated/cpp/type", "generated/cpp/nested/type", "quenching"):
                def pathto(docname):
                    if docname == pagename:
                        return "#"
                    return posixpath.relpath(docname + ".html", posixpath.dirname(pagename) or ".")
                rendered = template.render(cpp_navigation=navigation, pathto=pathto, pagename=pagename)
                self.assertIn("<h3>Python</h3>", rendered)
                self.assertIn("<h3>C++</h3>", rendered)
                for entry in navigation:
                    href = pathto(entry["docname"])
                    if entry["anchor"]:
                        href = href.rstrip("#") + "#" + entry["anchor"]
                    self.assertIn(f'href="{href}"', rendered)
                self.assertNotIn('href="##', rendered)
                for docname in ("installation", "getting_started", "quenching", "potentials", "system_class",
                                "landscape_module", "database", "utils", "generated/source_api/index"):
                    self.assertIn(f'href="{pathto(docname)}"', rendered)
                if pagename in ("generated/cpp/index", "quenching"):
                    self.assertEqual(rendered.count('aria-current="page"'), 1)

    def test_builder_hook_runs_after_exhale_and_compound_pages(self):
        extension = self.extension()
        connections = []
        extension.setup(SimpleNamespace(connect=lambda *args, **kwargs: connections.append((args, kwargs))))
        self.assertEqual(connections, [(("builder-inited", extension.generate), {"priority": 950})])


if __name__ == "__main__":
    unittest.main()

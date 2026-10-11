from pathlib import Path
import json
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

from docutils import nodes

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import conf


class DoxygenReferenceTests(unittest.TestCase):
    def test_missing_breathe_targets_link_to_full_doxygen_reference(self):
        with TemporaryDirectory() as directory:
            Path(directory, "index.xml").write_text('''<doxygenindex>
                <compound refid="classpele_1_1Array"><name>pele::Array</name></compound>
                <compound refid="namespacepele"><name>pele</name>
                    <member refid="namespacepele_1a123" kind="enum"><name>StopCriterionType</name></member>
                    <member refid="namespacepele_1a123a456" kind="enumvalue"><name>GRADIENT</name></member>
                </compound>
            </doxygenindex>''')
            callbacks = {}
            app = SimpleNamespace(
                connect=lambda event, callback, **kwargs: callbacks.update({callback.__name__: callback}),
                builder=SimpleNamespace(format="html", get_target_uri=lambda doc: f"{doc}.html"),
            )
            with patch.dict(conf.breathe_projects, {"pele": directory}):
                conf.setup(app)
            resolve = callbacks["doxygen_reference"]
            for domain, target, expected in [
                ("cpp", "Array<double>", "../doxygen/classpele_1_1Array.html"),
                ("std", "classpele_1_1Array", "../doxygen/classpele_1_1Array.html"),
                ("cpp", "StopCriterionType::GRADIENT", "../doxygen/namespacepele.html#a123a456"),
            ]:
                node = {"refdomain": domain, "reftype": "ref", "reftarget": target, "refdoc": "generated/example"}
                result = resolve(app, None, node, nodes.Text(target))
                self.assertEqual(result["refuri"], expected)
            node["reftarget"] = "UnknownType"
            self.assertIsNone(resolve(app, None, node, nodes.Text("UnknownType")))
            signature = callbacks["dict_method_signature"](
                app, "method", "Result.update", dict.update, None, None, None
            )
            self.assertEqual(signature, ("([E, ]**F)", "None"))
            for method in (dict.clear, dict.copy, dict.items, dict.keys, dict.values):
                signature = callbacks["dict_method_signature"](
                    app, "method", f"Result.{method.__name__}", method, None,
                    "()", "prose describing the result",
                )
                self.assertEqual(signature, ("()", None))

    def test_source_references_use_context_and_leave_ambiguous_names_unresolved(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "index.xml").write_text("<doxygenindex/>")
            output = root / "generated/source_api"
            output.mkdir(parents=True)
            entries = [
                {"fullname": "pele.example.Type.run", "docname": "api", "anchor": "method"},
                {"fullname": "pele.other.run", "docname": "other", "anchor": "function"},
            ]
            (output / "manifest.json").write_text(json.dumps(entries))
            callbacks = {}
            app = SimpleNamespace(
                srcdir=directory,
                connect=lambda event, callback, **kwargs: callbacks.update({callback.__name__: callback}),
                builder=SimpleNamespace(format="html", get_relative_uri=lambda source, target: target + ".html"),
            )
            with patch.dict(conf.breathe_projects, {"pele": directory}):
                conf.setup(app)
            resolve = callbacks["source_reference"]
            node = {"refdomain": "py", "reftarget": "run", "refdoc": "index",
                    "py:module": "pele.example", "py:class": "Type"}
            result = resolve(app, None, node, nodes.Text("run"))
            self.assertEqual(result["refuri"], "api.html#method")
            node.update({"py:module": None, "py:class": None})
            self.assertIsNone(resolve(app, None, node, nodes.Text("run")))
            node["reftarget"] = "pele.example.Type.run"
            result = resolve(app, None, node, nodes.Text("run"))
            self.assertEqual(result["refuri"], "api.html#method")

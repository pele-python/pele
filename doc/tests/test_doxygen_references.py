from pathlib import Path
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
                connect=lambda event, callback, **kwargs: callbacks.update({event: callback}),
                builder=SimpleNamespace(format="html", get_target_uri=lambda doc: f"{doc}.html"),
            )
            with patch.dict(conf.breathe_projects, {"pele": directory}):
                conf.setup(app)
            resolve = callbacks["missing-reference"]
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
            signature = callbacks["autodoc-process-signature"](
                app, "method", "Result.update", dict.update, None, None, None
            )
            self.assertEqual(signature, ("([E, ]**F)", "None"))

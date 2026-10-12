from pathlib import Path
import json
import posixpath
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

from docutils import nodes

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import conf


class DoxygenReferenceTests(unittest.TestCase):
    def test_doxygen_references_link_to_canonical_sphinx_cpp_pages(self):
        with TemporaryDirectory() as directory:
            Path(directory, "index.xml").write_text('''<doxygenindex>
                <compound refid="classpele_1_1Array"><name>pele::Array</name></compound>
                <compound refid="namespacepele"><name>pele</name>
                    <member refid="namespacepele_1a123" kind="enum"><name>StopCriterionType</name></member>
                    <member refid="namespacepele_1a123a456" kind="enumvalue"><name>GRADIENT</name></member>
                </compound>
                <compound refid="xsum_8h" kind="file"><name>xsum.h</name>
                    <member refid="xsum_8h_1typedef" kind="typedef"><name>xsum_small_accumulator</name></member>
                </compound>
                <compound refid="group__Enumerations" kind="group"><name>Enumerations</name>
                    <member refid="group__Enumerations_1ga123" kind="enum"><name>Termination</name></member>
                    <member refid="group__Enumerations_1gga123a456" kind="enumvalue"><name>ARMIJO</name></member>
                </compound>
                <compound refid="classpeleNode" kind="class"><name>pele::Node</name>
                    <member refid="classpeleNode_1typedef" kind="typedef"><name>edge_list</name></member>
                </compound>
                <compound refid="structpeleFactorial" kind="struct"><name>pele::Factorial</name>
                    <member refid="structpeleFactorial_1enum" kind="enum"><name /></member>
                    <member refid="structpeleFactorial_1value" kind="enumvalue"><name>value</name></member>
                </compound>
            </doxygenindex>''')
            for identifier, members in {
                "classpele_1_1Array": "",
                "namespacepele": '''<memberdef id="namespacepele_1a123" kind="enum">
                    <qualifiedname>pele::StopCriterionType</qualifiedname>
                    <enumvalue id="namespacepele_1a123a456"><name>GRADIENT</name></enumvalue></memberdef>''',
                "xsum_8h": "",
                "group__Enumerations": '''<memberdef id="group__Enumerations_1ga123" kind="enum">
                    <qualifiedname>pele::Termination</qualifiedname>
                    <enumvalue id="group__Enumerations_1gga123a456"><name>ARMIJO</name></enumvalue></memberdef>''',
                "classpeleNode": '''<memberdef id="classpeleNode_1typedef" kind="typedef">
                    <qualifiedname>pele::Node::edge_list</qualifiedname></memberdef>''',
                "structpeleFactorial": '''<memberdef id="structpeleFactorial_1enum" kind="enum">
                    <enumvalue id="structpeleFactorial_1value"><name>value</name></enumvalue></memberdef>''',
            }.items():
                Path(directory, identifier + ".xml").write_text(
                    f"<doxygen><compounddef><sectiondef>{members}</sectiondef></compounddef></doxygen>")
            callbacks = {}
            app = SimpleNamespace(
                connect=lambda event, callback, **kwargs: callbacks.update({callback.__name__: callback}),
                srcdir=directory,
                exhale_root=SimpleNamespace(root_directory=str(Path(directory, "generated/cpp")), node_by_refid={
                    "classpele_1_1Array": SimpleNamespace(file_name=str(Path(directory, "generated/cpp/array.rst")), parent=None),
                    "namespacepele": SimpleNamespace(file_name=str(Path(directory, "generated/cpp/namespace.rst")), parent=None),
                    "namespacepele_1a123": SimpleNamespace(file_name=str(Path(directory, "generated/cpp/enum.rst")), parent=None),
                    "xsum_8h": SimpleNamespace(file_name="xsum.rst", parent=None),
                    "group__Enumerations_1ga123": SimpleNamespace(file_name="termination.rst", parent=None),
                    "classpeleNode": SimpleNamespace(file_name="node.rst", parent=None),
                    "structpeleFactorial": SimpleNamespace(file_name="factorial.rst", parent=None),
                }),
                builder=SimpleNamespace(format="html", get_relative_uri=lambda source, target:
                    posixpath.relpath(target + ".html", posixpath.dirname(source) or ".")),
            )
            with patch.dict(conf.breathe_projects, {"pele": directory}):
                conf.setup(app)
            resolve = callbacks["doxygen_reference"]
            for domain, target, expected in [
                ("cpp", "Array<double>", "cpp/array.html#classpele_1_1Array"),
                ("std", "classpele_1_1Array", "cpp/array.html#classpele_1_1Array"),
                ("cpp", "StopCriterionType::GRADIENT", "cpp/enum.html#namespacepele_1a123a456"),
                ("cpp", "xsum_small_accumulator", "cpp/xsum.html#xsum_8h_1typedef"),
                ("c", "xsum_small_accumulator", "cpp/xsum.html#xsum_8h_1typedef"),
                ("cpp", "ARMIJO", "cpp/termination.html#group__Enumerations_1gga123a456"),
                ("std", "group__Enumerations_1gga123a456", "cpp/termination.html#group__Enumerations_1gga123a456"),
                ("cpp", "edge_list::iterator", "cpp/node.html#classpeleNode_1typedef"),
                ("cpp", "Factorial<m_p>::value", "cpp/factorial.html#structpeleFactorial_1value"),
            ]:
                node = {"refdomain": domain, "reftype": "ref", "reftarget": target, "refdoc": "generated/example"}
                result = resolve(app, None, node, nodes.Text(target))
                self.assertEqual(result["refuri"], expected)
            node["reftarget"] = "UnknownType"
            self.assertIsNone(resolve(app, None, node, nodes.Text("UnknownType")))
            app.builder.get_relative_uri = lambda source, target: posixpath.relpath(target, source) + "/"
            node["reftarget"] = "classpele_1_1Array"
            node["refdomain"] = "std"
            self.assertEqual(resolve(app, None, node, nodes.Text("Array"))["refuri"],
                             "../cpp/array/#classpele_1_1Array")
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

    def test_python312_inherited_enum_docstrings_render_byteorder_literal(self):
        from io import StringIO
        from docutils.core import publish_doctree
        with TemporaryDirectory() as directory:
            Path(directory, "index.xml").write_text("<doxygenindex/>")
            callbacks = {}
            app = SimpleNamespace(connect=lambda event, callback, **kwargs: callbacks.update({callback.__name__: callback}))
            with patch.dict(conf.breathe_projects, {"pele": directory}):
                conf.setup(app)
            callback = callbacks["int_bytes_docstring"]
            for enum in ("HessianType", "StopCriterionType"):
                for method in ("from_bytes", "to_bytes"):
                    lines = ["To request the native byte order, use `sys.byteorder' as the byte order value."]
                    callback(app, "method", f"pele.optimize.{enum}.{method}", None, None, lines)
                    tree = publish_doctree("\n".join(lines), settings_overrides={"warning_stream": StringIO()})
                    self.assertFalse(list(tree.findall(nodes.system_message)))
                    self.assertIn("``sys.byteorder``", lines[0])
            lines = ["Do not change another object's docstring."]
            callback(app, "method", "Example.to_bytes", None, None, lines)
            self.assertEqual(lines, ["Do not change another object's docstring."])

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

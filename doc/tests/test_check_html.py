from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_html import broken_links, check, coverage_errors, documented_cpp_symbols, sphinx_cpp_targets


class CheckHtmlTests(unittest.TestCase):
    def test_source_coverage_requires_the_rendered_anchor(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "api.html").write_text('<span id="documented"></span>')
            entries = [
                {"fullname": "covered", "docname": "api", "anchor": "documented"},
                {"fullname": "omitted", "docname": "api", "anchor": "missing"},
                {"fullname": "outside", "docname": "../api", "anchor": "documented"},
            ]
            self.assertEqual(coverage_errors(root, entries, "source API"), [
                "Missing source API target: omitted", "Missing source API page: outside",
            ])

    def test_native_coverage_keeps_developer_apis_and_excludes_private_helpers(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "api.xml").write_text('''<doxygen><compounddef id="classExample" kind="class">
                <compoundname>Example</compoundname><briefdescription> </briefdescription>
                <sectiondef>
                  <memberdef id="classExample_1aprivate" kind="function" prot="private">
                    <name>helper</name><briefdescription><para>Private helper.</para></briefdescription>
                  </memberdef>
                  <memberdef id="classExample_1adeveloper" kind="function" prot="protected">
                    <name>extension_point</name><briefdescription><para>Developer interface.</para></briefdescription>
                  </memberdef>
                  <memberdef id="classExample_1astatic" kind="function" prot="public" static="yes">
                    <name>factory</name><briefdescription><para>Public factory.</para></briefdescription>
                  </memberdef>
                  <memberdef id="file_1amacro" kind="define">
                    <name>FLAG</name><detaileddescription><para>File documentation.</para></detaileddescription>
                  </memberdef>
                  <memberdef id="file_1aempty" kind="function">
                    <name>undocumented</name><briefdescription> </briefdescription>
                  </memberdef>
                </sectiondef>
                </compounddef></doxygen>''')
            (root / "file.xml").write_text('''<doxygen><compounddef id="file" kind="file">
                <compoundname>implementation.cpp</compoundname><sectiondef>
                  <memberdef id="file_1alocal" kind="function" static="yes">
                    <name>local_helper</name><briefdescription><para>Private implementation.</para></briefdescription>
                  </memberdef>
                </sectiondef></compounddef></doxygen>''')
            (root / "stl.xml").write_text('''<doxygen><compounddef id="namespacestd" kind="namespace">
                <compoundname>std</compoundname><briefdescription><para>STL stub.</para></briefdescription>
                <location file="[STL]"/>
                </compounddef></doxygen>''')
            entries = documented_cpp_symbols(root)
            self.assertEqual([(item["fullname"], item["refid"]) for item in entries], [
                ("extension_point", "classExample_1adeveloper"),
                ("factory", "classExample_1astatic"),
                ("FLAG", "file_1amacro"),
            ])

    def test_native_targets_use_exact_generated_anchors_and_reject_missing_or_duplicate_pages(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            native = root / "generated/cpp"
            native.mkdir(parents=True)
            (native / "type.html").write_text('''<span id="classExample"></span>
                <span id="classExample_1arun"></span><span id="duplicate"></span>''')
            (native / "other.html").write_text('<span id="duplicate"></span>')
            (root / "cpp_reference.html").write_text('<span id="missing"></span>')
            entries = [
                {"fullname": "Example", "kind": "class", "refid": "classExample"},
                {"fullname": "Example::run", "kind": "function", "refid": "classExample_1arun"},
                {"fullname": "omitted", "kind": "function", "refid": "missing"},
                {"fullname": "repeated", "kind": "function", "refid": "duplicate"},
            ]
            rendered, errors = sphinx_cpp_targets(root, entries)
            self.assertEqual([(entry["refid"], entry["docname"], entry["anchor"]) for entry in rendered], [
                ("classExample", "generated/cpp/type", "classExample"),
                ("classExample_1arun", "generated/cpp/type", "classExample_1arun"),
            ])
            self.assertEqual(errors, [
                "Missing Sphinx C++ API target: omitted (missing)",
                "Duplicate Sphinx C++ API target: repeated (duplicate): generated/cpp/other.html, generated/cpp/type.html",
            ])
            self.assertEqual(coverage_errors(root, rendered, "Sphinx C++ API"), [])

    def test_links_and_anchors(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "index.html").write_text('''
                <a href="api.html#method">API</a>
                <a href="api.html#missing">Broken anchor</a>
                <a href="absent.html">Broken page</a>
                <a href="../index.html">Outside the site</a>
                <img src="missing.png">
                <a href="https://example.com/">External</a>
                <a href="//example.com/">External</a>
                <a href="mailto:user@example.com">Mail</a>
            ''')
            (root / "api.html").write_text('<a name="method"></a><a href="/index.html">Home</a>')
            self.assertEqual(broken_links(root), [
                "index.html: missing anchor api.html#missing",
                "index.html: missing file absent.html",
                "index.html: missing file ../index.html",
                "index.html: missing file missing.png",
            ])

    def test_missing_generation_fails(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            with patch("check_html.documented_cpp_symbols", return_value=[]):
                with self.assertRaises(SystemExit) as result:
                    check(root)
            self.assertIn("Missing generated page: index.html", str(result.exception))
            self.assertIn("Missing generated page: generated/cpp/index.html", str(result.exception))


if __name__ == "__main__":
    unittest.main()

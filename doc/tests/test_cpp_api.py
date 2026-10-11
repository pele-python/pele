import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import xml.etree.ElementTree as ET


def extension():
    source = Path(__file__).resolve().parents[1] / "sphinxext/cpp_api.py"
    spec = importlib.util.spec_from_file_location("cpp_api", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CppApiTests(unittest.TestCase):
    def test_global_declarations_share_one_page_and_keep_complementary_documentation(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            functions = {
                "header": ("source/pele/example.h", "source/example.c", "double", "header_function", "Declaration docs."),
                "implementation": ("source/example.c", "source/example.c", "double", "body_function", "Implementation docs."),
                "different_body": ("source/other.c", "source/other.c", "double", "other_function", "Other body docs."),
                "overload": ("source/pele/example.h", "source/example.c", "int", "overload_function", "Overload docs."),
            }
            index = []
            for identifier, (filename, bodyfile, parameter_type, member_id, description) in functions.items():
                index.append(f'<compound kind="file" refid="{identifier}"><name>{filename}</name>'
                             f'<member kind="function" refid="{member_id}"><name>example</name></member></compound>')
                declaration = f'declfile="{filename}"' if identifier == "header" else ""
                extra = "<para>Body details.</para>" if identifier == "implementation" else ""
                (root / f"{identifier}.xml").write_text(f'''<doxygen><compounddef id="{identifier}" kind="file">
                    <compoundname>{filename}</compoundname><sectiondef><memberdef id="{member_id}" kind="function" static="no" prot="public">
                    <type>void</type><name>example</name><definition>void example</definition>
                    <param><type>{parameter_type}</type><declname>{'' if identifier == 'header' else 'argument'}</declname></param>
                    <briefdescription><para>{description}</para></briefdescription>
                    <detaileddescription><para>Shared docs.</para>{extra}</detaileddescription>
                    <location file="{filename}" bodyfile="{bodyfile}" bodystart="10" bodyend="20" {declaration}/>
                    </memberdef></sectiondef><listofallmembers><member refid="{member_id}"/></listofallmembers>
                    <programlisting><codeline><ref refid="body_function">example</ref></codeline></programlisting>
                    </compounddef></doxygen>''')
            (root / "index.xml").write_text("<doxygenindex>" + "".join(index) + "</doxygenindex>")
            module = extension()
            module.prune_xml(root)
            header = ET.parse(root / "header.xml").getroot()
            member = header.find(".//memberdef")
            descriptions = " ".join(member.itertext())
            for text in ("Declaration docs.", "Implementation docs.", "Shared docs.", "Body details."):
                self.assertIn(text, descriptions)
            self.assertEqual(descriptions.count("Shared docs."), 1)
            self.assertEqual(ET.parse(root / "implementation.xml").getroot().findall(".//memberdef"), [])
            self.assertEqual(ET.parse(root / "implementation.xml").getroot().findall(".//listofallmembers/member"), [])
            for identifier in ("different_body", "overload"):
                self.assertIsNotNone(ET.parse(root / f"{identifier}.xml").getroot().find(".//memberdef"))
            index_members = [node.get("refid") for node in ET.parse(root / "index.xml").getroot().findall(".//member")]
            self.assertEqual(index_members, ["header_function", "other_function", "overload_function"])
            for path in root.glob("*.xml"):
                self.assertNotIn('refid="body_function"', path.read_text())
            before = {path: path.read_bytes() for path in root.glob("*.xml")}
            module.prune_xml(root)
            self.assertEqual(before, {path: path.read_bytes() for path in root.glob("*.xml")})

    def test_generated_paths_use_the_same_project_relative_directory_and_file_names(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            project = Path("/checkout/pele")
            (root / "index.xml").write_text('''<doxygenindex>
                <compound kind="dir" refid="dirSource"><name>/checkout/pele/source</name></compound>
                <compound kind="file" refid="fileExample"><name>/checkout/pele/source/pele/example.h</name></compound>
                </doxygenindex>''')
            (root / "dirSource.xml").write_text('''<doxygen><compounddef id="dirSource" kind="dir">
                <compoundname>/checkout/pele/source</compoundname>
                <innerdir refid="dirPele">/checkout/pele/source/pele</innerdir>
                <innerfile refid="fileExample">/checkout/pele/source/pele/example.h</innerfile>
                <location file="/checkout/pele/source"/></compounddef></doxygen>''')
            (root / "fileExample.xml").write_text('''<doxygen><compounddef id="fileExample" kind="file">
                <compoundname>example.h</compoundname>
                <location file="/checkout/pele/source/pele/example.h"/>
                <sectiondef><memberdef id="example" kind="function"><name>example</name>
                <location file="/checkout/pele/source/pele/example.h" bodyfile="/checkout/pele/source/example.c"
                  declfile="/checkout/pele/source/pele/example.h"/></memberdef></sectiondef>
                </compounddef></doxygen>''')
            module = extension()
            module.prune_xml(root, project_root=project)
            for path in root.glob("*.xml"):
                self.assertNotIn(str(project), path.read_text())
            directory = ET.parse(root / "dirSource.xml").getroot().find("compounddef")
            self.assertEqual(directory.findtext("compoundname"), "source")
            self.assertEqual(directory.findtext("innerdir"), "source/pele")
            location = ET.parse(root / "fileExample.xml").getroot().find(".//memberdef/location")
            self.assertEqual(location.attrib, {"file": "source/pele/example.h", "bodyfile": "source/example.c", "declfile": "source/pele/example.h"})
            self.assertRegex((root / "fileExample.xml").read_text(),
                             r'<location file="source/pele/example.h"/>')

    def test_only_identical_qualified_aliases_are_deduplicated_and_references_remapped(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            names = {
                "classArray": ("Array", "array.hpp", "pele::Array<T>::size"),
                "classpeleArray": ("pele::Array", "array.hpp", "pele::Array<T>::size"),
                "classOtherArray": ("other::Array", "other.hpp", "pele::Array<T>::size"),
                "classSpecialized": ("pele::Array<int>", "array.hpp", "pele::Array<T>::size"),
                "classDifferent": ("different::Array", "array.hpp", "different::Array<T>::size"),
            }
            compounds = []
            for identifier, (name, filename, definition) in names.items():
                compounds.append(f'<compound kind="class" refid="{identifier}"><name>{name}</name></compound>')
                (root / f"{identifier}.xml").write_text(f'''<doxygen><compounddef id="{identifier}" kind="class" prot="public">
                    <compoundname>{name.replace('<', '&lt;').replace('>', '&gt;')}</compoundname>
                    <location file="{filename}" line="10" column="1" bodyfile="{filename}" bodystart="10" bodyend="20"/>
                    <sectiondef><memberdef id="{identifier}_1size" kind="function" prot="public">
                    <name>size</name><definition>{definition.replace('<', '&lt;').replace('>', '&gt;')}</definition>
                    <argsstring>() const</argsstring></memberdef></sectiondef>
                    </compounddef></doxygen>''')
            (root / "index.xml").write_text("<doxygenindex>" + "".join(compounds).replace("<int>", "&lt;int&gt;") + "</doxygenindex>")
            (root / "file.xml").write_text('''<doxygen><compounddef id="file" kind="file"><compoundname>example.cpp</compoundname>
                <innerclass refid="classArray">Array</innerclass>
                <briefdescription><para><ref refid="classArray_1size">Array::size</ref></para></briefdescription>
                <programlisting><codeline><highlight><ref refid="classArray">Array</ref></highlight></codeline></programlisting>
                </compounddef></doxygen>''')
            module = extension()
            module.prune_xml(root)
            self.assertFalse((root / "classArray.xml").exists())
            for identifier in names.keys() - {"classArray"}:
                self.assertTrue((root / f"{identifier}.xml").exists())
            index_ids = [node.get("refid") for node in ET.parse(root / "index.xml").getroot()]
            self.assertEqual(index_ids, [identifier for identifier in names if identifier != "classArray"])
            references = [node.get("refid") for node in ET.parse(root / "file.xml").getroot().iter() if node.get("refid")]
            self.assertEqual(references, ["classpeleArray_1size", "classpeleArray"])
            self.assertEqual(ET.parse(root / "file.xml").getroot().findall(".//innerclass"), [])
            before = {path: path.read_bytes() for path in root.glob("*.xml")}
            module.prune_xml(root)
            self.assertEqual(before, {path: path.read_bytes() for path in root.glob("*.xml")})

    def test_generated_xml_omits_helpers_and_preserves_developer_interfaces(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "index.xml").write_text('''<doxygenindex>
                <compound kind="class" refid="classExample"><name>Example</name>
                  <member kind="function" refid="private"><name>hidden</name></member>
                  <member kind="function" refid="factory"><name>factory</name></member>
                  <member kind="function" refid="protected"><name>extension_point</name></member>
                </compound>
                <compound kind="namespace" refid="namespacepele"><name>pele</name>
                  <member kind="function" refid="local"><name>local_helper</name></member>
                  <member kind="function" refid="global"><name>public_function</name></member>
                </compound>
                <compound kind="class" refid="classPrivate"><name>_Private</name></compound>
                <compound kind="namespace" refid="namespacestd"><name>std</name></compound>
                </doxygenindex>''')
            (root / "classExample.xml").write_text('''<doxygen><compounddef id="classExample" kind="class" prot="public">
                <compoundname>Example</compoundname><innerclass refid="classPrivate">_Private</innerclass>
                <sectiondef><memberdef id="private" kind="function" prot="private"><name>hidden</name></memberdef>
                  <memberdef id="underscore" kind="function" prot="public"><name>_helper</name></memberdef>
                  <memberdef id="factory" kind="function" prot="public" static="yes"><name>factory</name></memberdef>
                  <memberdef id="protected" kind="function" prot="protected"><name>extension_point</name></memberdef>
                  <memberdef id="privateEnum" kind="enum" prot="private"><name>PrivateEnum</name>
                    <enumvalue id="privateValue"><name>VALUE</name></enumvalue></memberdef>
                </sectiondef><briefdescription><para>Before <ref refid="private">hidden helper</ref> after.</para></briefdescription>
                <listofallmembers><member refid="private"/></listofallmembers>
                </compounddef></doxygen>''')
            (root / "namespacepele.xml").write_text('''<doxygen><compounddef id="namespacepele" kind="namespace">
                <compoundname>pele</compoundname><innerclass refid="classPrivate">_Private</innerclass>
                <sectiondef><memberdef id="local" kind="function" static="yes"><name>local_helper</name></memberdef>
                  <memberdef id="global" kind="function" static="no"><name>public_function</name></memberdef>
                </sectiondef></compounddef></doxygen>''')
            (root / "classPrivate.xml").write_text('''<doxygen><compounddef id="classPrivate" kind="class">
                <compoundname>_Private</compoundname><sectiondef>
                <memberdef id="privateOwnerMember" kind="function"><name>public_spelling</name></memberdef>
                </sectiondef></compounddef></doxygen>''')
            (root / "namespacestd.xml").write_text('''<doxygen><compounddef id="namespacestd" kind="namespace">
                <compoundname>std</compoundname><location file="[STL]"/></compounddef></doxygen>''')
            module = extension()
            excluded = module.prune_xml(root)
            self.assertEqual(excluded, {"private", "underscore", "privateEnum", "privateValue", "local",
                                        "classPrivate", "privateOwnerMember", "namespacestd"})
            index = ET.parse(root / "index.xml").getroot()
            self.assertEqual([node.get("refid") for node in index.iter() if node.get("refid")],
                             ["classExample", "factory", "protected", "namespacepele", "global"])
            example = ET.parse(root / "classExample.xml").getroot()
            self.assertEqual([node.get("id") for node in example.findall(".//memberdef")], ["factory", "protected"])
            self.assertEqual(example.findtext(".//para"), "Before hidden helper after.")
            self.assertEqual(example.findall(".//innerclass"), [])
            self.assertFalse((root / "classPrivate.xml").exists())
            self.assertFalse((root / "namespacestd.xml").exists())
            self.assertEqual(module.prune_xml(root), set())

    def test_aliases_keep_exact_ids_and_group_descriptions_do_not_repeat_declarations(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            xml = root / "xml"
            output = root / "generated/cpp"
            xml.mkdir()
            output.mkdir(parents=True)
            nodes = {}
            for identifier, kind in [("File_MixedCase", "file"), ("namespacepele", "namespace"), ("group__Enumerations", "group")]:
                (xml / f"{identifier}.xml").write_text(f'''<doxygen><compounddef id="{identifier}" kind="{kind}">
                    <compoundname>{identifier}</compoundname><title>{identifier}</title>
                    <briefdescription><para>A public description.</para></briefdescription>
                    <sectiondef><memberdef id="enum" kind="enum"><name>Flags</name></memberdef></sectiondef>
                    </compounddef></doxygen>''')
                nodes[identifier] = SimpleNamespace(file_name=f"{identifier}.rst", link_name=f"exhale_{identifier}")
                if kind != "group":
                    (output / f"{identifier}.rst").write_text(f"{identifier}\n================\n\nExisting description.\n")
            (output / "index.rst").write_text("C++ API\n=======\n")
            exhale_root = SimpleNamespace(root_directory=str(output), full_root_file_path=str(output / "index.rst"), node_by_refid=nodes)
            module = extension()
            aliases = module.add_compound_pages(xml, exhale_root)
            self.assertEqual(set(aliases), set(nodes))
            before = {path: path.read_text() for path in output.glob("*.rst")}
            module.add_compound_pages(xml, exhale_root)
            self.assertEqual(before, {path: path.read_text() for path in output.glob("*.rst")})
            for identifier in nodes:
                contents = (output / f"{identifier}.rst").read_text()
                self.assertEqual(contents.count(f'<span id="{identifier}">'), 1)
                self.assertNotIn(f".. _{identifier}:", contents)
            group = (output / "group__Enumerations.rst").read_text()
            self.assertIn(".. doxygengroup:: group__Enumerations", group)
            self.assertIn(":desc-only:", group)
            self.assertIn(":no-link:", group)
            self.assertNotIn(".. doxygenenum::", group)
            self.assertNotIn(":members:", group)
            self.assertIn("   group__Enumerations\n", (output / "index.rst").read_text())


if __name__ == "__main__":
    unittest.main()

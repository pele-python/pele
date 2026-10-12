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
    def test_global_types_variables_and_macros_keep_one_canonical_documented_definition(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            index = []
            for identifier, filename, variable_initializer in (
                    ("header", "source/pele/example.hpp", ""),
                    ("second", "source/pele/second.hpp", ""),
                    ("body", "source/example.cpp", "= 0")):
                records = f'''<memberdef id="{identifier}_Scalar" kind="typedef"><type>double</type><name>Scalar</name>
                    <definition>using Scalar = double</definition><briefdescription><para>{identifier} type docs.</para></briefdescription>
                    <location file="{filename}"/></memberdef>
                    <memberdef id="{identifier}_debug" kind="variable"><type>int</type><name>debug</name>
                    <definition>int debug</definition><initializer>{variable_initializer}</initializer>
                    <briefdescription><para>{identifier} variable docs.</para></briefdescription><location file="{filename}"/></memberdef>
                    <memberdef id="{identifier}_mode" kind="define"><name>MODE</name><initializer>1</initializer>
                    <briefdescription><para>{identifier} macro docs.</para></briefdescription><location file="{filename}"/></memberdef>'''
                if identifier == "header":
                    records += '''<memberdef id="restrict" kind="define"><name>restrict</name><initializer>__restrict__</initializer></memberdef>
                        <memberdef id="public_macro" kind="define"><name>PUBLIC_MACRO</name><briefdescription><para>Public API.</para></briefdescription></memberdef>
                        <memberdef id="different_type" kind="typedef"><type>int</type><name>Scalar</name><definition>using Scalar = int</definition></memberdef>'''
                (root / f"{identifier}.xml").write_text(f'''<doxygen><compounddef id="{identifier}" kind="file">
                    <compoundname>{filename}</compoundname><sectiondef>{records}</sectiondef></compounddef></doxygen>''')
                members = ''.join(f'<member kind="{kind}" refid="{identifier}_{name}"><name>{name}</name></member>'
                                  for kind, name in (("typedef", "Scalar"), ("variable", "debug"), ("define", "mode")))
                index.append(f'<compound refid="{identifier}" kind="file"><name>{filename}</name>{members}'
                             '<member kind="function" refid="stale"><name>missing_static_helper</name></member></compound>')
            (root / "index.xml").write_text('<doxygenindex>' + ''.join(index) + '</doxygenindex>')
            module = extension()
            module.prune_xml(root)
            members = [member for path in root.glob('*.xml') for member in ET.parse(path).getroot().findall('.//memberdef')]
            self.assertEqual(sum(member.findtext('definition') == 'using Scalar = double' for member in members), 1)
            variable = [member for member in members if member.findtext('name') == 'debug']
            self.assertEqual(len(variable), 1)
            self.assertEqual(variable[0].findtext('initializer'), '= 0')
            for identifier in ('header', 'second', 'body'):
                self.assertIn(f'{identifier} variable docs.', ''.join(variable[0].itertext()))
            self.assertEqual(sum(member.findtext('name') == 'MODE' for member in members), 1)
            self.assertIn('public_macro', {member.get('id') for member in members})
            self.assertIn('different_type', {member.get('id') for member in members})
            self.assertNotIn('restrict', {member.get('id') for member in members})
            self.assertNotIn('refid="stale"', (root / 'index.xml').read_text())

    def test_anonymous_enums_keep_values_and_drop_only_empty_foreign_copies(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for identifier, copied in (("structPrimary", ""), ("structSpecialized", '<memberdef kind="enum" id="structPrimary_1enum"><name/></memberdef>')):
                (root / f'{identifier}.xml').write_text(f'''<doxygen><compounddef id="{identifier}" kind="struct">
                    <compoundname>{identifier}</compoundname><sectiondef>{copied}
                    <memberdef kind="enum" id="{identifier}_1enum"><name/><enumvalue id="{identifier}_1value"><name>value</name></enumvalue></memberdef>
                    <memberdef kind="enum" id="{identifier}_1empty"><name/></memberdef>
                    </sectiondef></compounddef></doxygen>''')
            module = extension()
            module.prune_xml(root)
            enums = [member for path in root.glob('*.xml') for member in ET.parse(path).getroot().findall('.//memberdef')]
            self.assertEqual(len(enums), 4)
            self.assertEqual(sum(member.get('id') == 'structPrimary_1enum' for member in enums), 1)
            self.assertEqual(sum(len(member.findall('enumvalue')) for member in enums), 2)
            for member in enums:
                self.assertEqual(member.findtext('name'), '@anonymous_' + member.get('id'))

    def test_generated_item_directives_use_exact_ids_and_heading_names_stay_literal(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            xml = root / 'xml'
            output = root / 'generated/cpp'
            xml.mkdir()
            output.mkdir(parents=True)
            (xml / 'namespacepele.xml').write_text('''<doxygen><compounddef id="namespacepele" kind="namespace"><compoundname>pele</compoundname>
                <sectiondef><memberdef id="group_Flags_1enum" kind="enum"><name>Flags</name></memberdef></sectiondef></compounddef></doxygen>''')
            (xml / 'classExample.xml').write_text('''<doxygen><compounddef id="classExample" kind="class"><compoundname>Example</compoundname>
                <sectiondef><memberdef id="classExample_1enum" kind="enum"><name>Inner</name></memberdef></sectiondef></compounddef></doxygen>''')
            (xml / 'groupFlags.xml').write_text('''<doxygen><compounddef id="groupFlags" kind="group"><compoundname>Flags</compoundname>
                <sectiondef><memberdef id="groupOnly_1enum" kind="enum"><name>GroupOnly</name><qualifiedname>pele::GroupOnly</qualifiedname></memberdef></sectiondef></compounddef></doxygen>''')
            (output / 'enum.rst').write_text('''Enum Flags\n==========\n\n.. doxygenenum:: pele::Flags\n   :project: pele\n''')
            (output / 'struct.rst').write_text('''Struct UserData_\n================\n\n.. doxygenstruct:: pele::UserData_\n''')
            (output / 'group-enum.rst').write_text('''Enum GroupOnly\n==============\n\n.. doxygenenum:: pele::GroupOnly\n''')
            exhale_root = SimpleNamespace(root_directory=str(output), node_by_refid={
                'group_Flags_1enum': SimpleNamespace(file_name='enum.rst', kind='enum'),
                'groupOnly_1enum': SimpleNamespace(file_name='group-enum.rst', kind='enum'),
                'classExample_1enum': SimpleNamespace(file_name=None, kind='enum')})
            module = extension()
            module.rewrite_item_pages(xml, exhale_root)
            self.assertIn('.. doxygenrefid:: group_Flags_1enum\n   :owner: namespacepele', (output / 'enum.rst').read_text())
            self.assertIn(':project: pele', (output / 'enum.rst').read_text())
            self.assertIn(':qualified-name: pele::GroupOnly', (output / 'group-enum.rst').read_text())
            title = (output / 'struct.rst').read_text().splitlines()
            self.assertEqual(title[0], r'Struct UserData\_')
            self.assertEqual(len(title[0]), len(title[1]))
            self.assertIn('.. doxygenstruct:: pele::UserData_', (output / 'struct.rst').read_text())
            before = {path: path.read_bytes() for path in output.glob('*.rst')}
            module.rewrite_item_pages(xml, exhale_root)
            self.assertEqual(before, {path: path.read_bytes() for path in output.glob('*.rst')})

    def test_global_declarations_share_one_page_and_keep_complementary_documentation(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            functions = {
                "header": ("source/pele/example.h", "source/example.c", "double", "header_function", "Declaration docs."),
                "implementation": ("source/example.c", "source/example.c", "double", "body_function", "Implementation docs."),
                "different_body": ("source/other.c", "source/other.c", "double", "other_function", "Other body docs."),
                "overload": ("source/pele/example.h", "source/example.c", "int", "overload_function", "Overload docs."),
                "restrict_header": ("source/pele/restrict.h", "source/restrict.c", "double *", "restrict_header_function", "Restricted declaration docs."),
                "restrict_body": ("source/restrict.c", "source/restrict.c", "double *restrict", "restrict_body_function", "Restricted implementation docs."),
                "ambiguous_header": ("source/pele/ambiguous.h", "", "float *", "ambiguous_header_function", "Ambiguous declaration docs."),
                "ambiguous_first": ("source/first.c", "source/first.c", "float *restrict", "ambiguous_first_function", "First distinct body."),
                "ambiguous_second": ("source/second.c", "source/second.c", "float *", "ambiguous_second_function", "Second distinct body."),
            }
            index = []
            for identifier, (filename, bodyfile, parameter_type, member_id, description) in functions.items():
                index.append(f'<compound kind="file" refid="{identifier}"><name>{filename}</name>'
                             f'<member kind="function" refid="{member_id}"><name>example</name></member></compound>')
                declaration = f'declfile="{filename}"' if identifier == "header" else ""
                extra = "<para>Body details.</para>" if identifier == "implementation" else ""
                body_location = '' if identifier in ('restrict_header', 'ambiguous_header') else f'bodyfile="{bodyfile}" bodystart="10" bodyend="20"'
                includes = '<includes refid="restrict_header">restrict.h</includes>' if identifier == 'restrict_body' else (
                    '<includes refid="ambiguous_header">ambiguous.h</includes>' if identifier.startswith('ambiguous_') and identifier != 'ambiguous_header' else '')
                (root / f"{identifier}.xml").write_text(f'''<doxygen><compounddef id="{identifier}" kind="file">
                    <compoundname>{filename}</compoundname>{includes}<sectiondef><memberdef id="{member_id}" kind="function" static="no" prot="public">
                    <type>void</type><name>example</name><definition>void example</definition>
                    <param><type>{parameter_type}</type><declname>{'' if identifier == 'header' else 'argument'}</declname></param>
                    <briefdescription><para>{description}</para></briefdescription>
                    <detaileddescription><para>Shared docs.</para>{extra}</detaileddescription>
                    <location file="{filename}" {body_location} {declaration}/>
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
            self.assertEqual(index_members, ["header_function", "other_function", "overload_function", "restrict_header_function",
                                             "ambiguous_header_function", "ambiguous_first_function", "ambiguous_second_function"])
            restricted = ET.parse(root / "restrict_header.xml").getroot().find(".//memberdef")
            self.assertEqual(restricted.findtext("param/type"), "double *")
            self.assertIn("Restricted implementation docs.", "".join(restricted.itertext()))
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


class ExactRefidRenderingTests(unittest.TestCase):
    def test_exact_refids_render_c_header_function_group_enum_and_macro(self):
        try:
            from sphinx.application import Sphinx
        except ImportError:
            self.skipTest("Sphinx/Breathe are needed for the rendering fixture")
        from io import StringIO
        CPP_API = Path(__file__).resolve().parents[1] / "sphinxext"
        with TemporaryDirectory(prefix='pele-refid-fixture-') as directory:
            root = Path(directory)
            source = root / 'src'
            xml = source / 'xml'
            xml.mkdir(parents=True)
            (source / 'conf.py').write_text(f'''
import sys
sys.path.insert(0, {str(CPP_API)!r})
extensions = ['breathe']
project = 'Exact native API IDs'
master_doc = 'index'
breathe_projects = {{'fixture': 'xml'}}
breathe_default_project = 'fixture'
breathe_domain_by_extension = {{'h': 'c', 'hpp': 'cpp'}}
def setup(app):
    from cpp_api import _refid_directive
    app.add_directive('doxygenrefid', _refid_directive())
''')
            (source / 'index.rst').write_text('''Exact native API IDs
====================

.. doxygenrefid:: file__fixture_8h_1c_restrict
   :owner: file__fixture_8h

.. doxygenrefid:: group__Flags_1Enum_MixedCase
   :owner: namespacepele

.. doxygenrefid:: file__fixture_8h_1PUBLIC_LIMIT
   :owner: file__fixture_8h

.. doxygenrefid:: group__Flags_1Enum_GroupOnly
   :owner: group__Flags
   :qualified-name: pele::GroupOnly
''')
            (xml / 'index.xml').write_text('''<doxygenindex version="1.9.1">
<compound kind="file" refid="file__fixture_8h"><name>fixture.h</name>
  <member kind="function" refid="file__fixture_8h_1c_restrict"><name>sum_values</name></member>
  <member kind="define" refid="file__fixture_8h_1PUBLIC_LIMIT"><name>PUBLIC_LIMIT</name></member>
</compound>
<compound kind="namespace" refid="namespacepele"><name>pele</name></compound>
<compound kind="group" refid="group__Flags"><name>Flags</name>
  <member kind="enum" refid="group__Flags_1Enum_MixedCase"><name>Flags</name></member>
</compound>
</doxygenindex>''')
            (xml / 'file__fixture_8h.xml').write_text('''<doxygen version="1.9.1"><compounddef id="file__fixture_8h" kind="file">
<compoundname>fixture.h</compoundname><location file="fixture.h"/>
<sectiondef kind="func"><memberdef kind="function" id="file__fixture_8h_1c_restrict" prot="public" static="no" const="no" explicit="no" inline="no" virt="non-virtual">
  <type>double</type><definition>double sum_values</definition><argsstring>(double *restrict values)</argsstring><name>sum_values</name>
  <param><type>double *restrict</type><declname>values</declname></param>
  <briefdescription><para>A documented C restrict function.</para></briefdescription>
  <detaileddescription><para>Sum the caller-owned values.</para></detaileddescription>
  <location file="fixture.h" line="1"/>
</memberdef></sectiondef>
<sectiondef kind="define"><memberdef kind="define" id="file__fixture_8h_1PUBLIC_LIMIT" prot="public">
  <name>PUBLIC_LIMIT</name><initializer>32</initializer>
  <briefdescription><para>The documented public header macro.</para></briefdescription><detaileddescription/>
  <location file="fixture.h" line="2"/>
</memberdef></sectiondef>
</compounddef></doxygen>''')
            enum = '''<memberdef kind="enum" id="group__Flags_1Enum_MixedCase" prot="public" strong="no">
  <type>int</type><name>Flags</name><qualifiedname>pele::Flags</qualifiedname>
  <briefdescription><para>The documented group-owned enum.</para></briefdescription><detaileddescription/>
  <enumvalue id="group__Flags_1Enum_MixedCase_1Public_Value" prot="public">
    <name>PUBLIC_VALUE</name><initializer>= 1</initializer>
    <briefdescription><para>A documented enum value.</para></briefdescription><detaileddescription/>
  </enumvalue>
  <location file="fixture.hpp" line="3"/>
</memberdef>'''
            (xml / 'namespacepele.xml').write_text('<doxygen version="1.9.1"><compounddef id="namespacepele" kind="namespace"><compoundname>pele</compoundname><sectiondef kind="enum">' + enum + '</sectiondef></compounddef></doxygen>')
            group_only = '''<memberdef kind="enum" id="group__Flags_1Enum_GroupOnly" prot="public" strong="no">
              <type>int</type><name>GroupOnly</name><qualifiedname>pele::GroupOnly</qualifiedname>
              <briefdescription><para>The enum exists only in the group.</para></briefdescription><detaileddescription/>
              <enumvalue id="group__Flags_1Enum_GroupOnly_1Ready" prot="public"><name>READY</name><initializer>= 1</initializer>
                <briefdescription><para>The group-only value.</para></briefdescription><detaileddescription/></enumvalue>
              <location file="fixture.hpp" line="4"/></memberdef>'''
            (xml / 'group__Flags.xml').write_text('<doxygen version="1.9.1"><compounddef id="group__Flags" kind="group"><compoundname>Flags</compoundname><title>Flags</title><sectiondef kind="enum">' + enum + group_only + '</sectiondef></compounddef></doxygen>')
            warnings = StringIO()
            app = Sphinx(str(source), str(source), str(root / 'html'), str(root / 'doctrees'),
                         'html', status=StringIO(), warning=warnings, freshenv=True,
                         warningiserror=True)
            app.build(force_all=True)
            self.assertEqual(warnings.getvalue(), '')
            self.assertEqual(app.statuscode, 0)
            html = (root / 'html/index.html').read_text()
            for identifier in ('file__fixture_8h_1c_restrict', 'group__Flags_1Enum_MixedCase',
                               'group__Flags_1Enum_MixedCase_1Public_Value', 'file__fixture_8h_1PUBLIC_LIMIT',
                               'group__Flags_1Enum_GroupOnly', 'group__Flags_1Enum_GroupOnly_1Ready'):
                self.assertEqual(html.count(f'id="{identifier}"'), 1, identifier)
            for description in ('A documented C restrict function.', 'Sum the caller-owned values.',
                                'The documented group-owned enum.', 'A documented enum value.',
                                'The documented public header macro.'):
                self.assertIn(description, html)
            self.assertIn('restrict', html)
            cpp_objects = {name for name, *_ in app.env.get_domain('cpp').get_objects()}
            self.assertIn('pele::GroupOnly', cpp_objects)
            self.assertIn('pele::GroupOnly::READY', cpp_objects)


if __name__ == "__main__":
    unittest.main()

"""Keep the complete Exhale reference public and provide exact Doxygen anchors."""

from copy import deepcopy
from html import escape
from pathlib import Path
import xml.etree.ElementTree as ET


def _private_name(name):
    return any(part.strip().startswith("_") for part in name.split("::"))


def _normalize_paths(trees, project_root):
    project_root = Path(project_root).resolve()

    def relative(value):
        path = Path(value)
        if path.is_absolute() and path.is_relative_to(project_root):
            return path.relative_to(project_root).as_posix()
        return value

    changed = set()
    for path, tree in trees.items():
        for node in tree.getroot().iter():
            if node.tag == "location":
                for attribute in ("file", "bodyfile", "declfile"):
                    value = node.get(attribute)
                    if value and relative(value) != value:
                        node.set(attribute, relative(value))
                        changed.add(path)
            elif node.tag in ("compoundname", "innerdir", "innerfile", "name") and node.text:
                value = node.text.strip()
                if relative(value) != value:
                    node.text = relative(value)
                    changed.add(path)
    return changed


def _text(node):
    return " ".join("".join(node.itertext()).split()) if node is not None else ""


def _description_key(node):
    return (node.tag, tuple(sorted(node.attrib.items())), " ".join((node.text or "").split()),
            tuple(_description_key(child) for child in node))


def _merge_descriptions(destination, source):
    for tag in ("briefdescription", "detaileddescription", "inbodydescription"):
        description = source.find(tag)
        if description is None or not _text(description):
            continue
        existing = destination.find(tag)
        if existing is None:
            destination.append(deepcopy(description))
            continue
        keys = {_description_key(child) for child in existing}
        for child in description:
            key = _description_key(child)
            if key not in keys:
                existing.append(deepcopy(child))
                keys.add(key)


def _deduplicate_functions(trees, excluded_files):
    """Join a global declaration and definition only when they share one body."""
    groups = {}
    for path, tree in trees.items():
        if path in excluded_files:
            continue
        compound = tree.getroot().find("compounddef")
        if compound is None or compound.get("kind") not in ("file", "namespace"):
            continue
        for member in compound.findall(".//memberdef"):
            location = member.find("location")
            if (member.get("kind") != "function" or member.get("static") == "yes"
                    or member.get("prot") == "private" or location is None
                    or not all(location.get(field) for field in ("bodyfile", "bodystart", "bodyend"))):
                continue
            signature = (
                member.findtext("qualifiedname") or member.findtext("name"), _text(member.find("type")),
                tuple(_text(param.find("type")) for param in member.findall("param")),
                tuple((_text(param.find("type")), _text(param.find("defval")))
                      for param in member.findall("templateparamlist/param")),
                tuple(location.get(field) for field in ("bodyfile", "bodystart", "bodyend")),
            )
            groups.setdefault(signature, []).append((path, member))
    aliases = {}
    changed = set()
    for records in groups.values():
        if len({member.get("id") for _, member in records}) < 2:
            continue

        def preference(record):
            location = record[1].find("location")
            header = Path(location.get("file", "")).suffix.lower() in (".h", ".hh", ".hpp", ".hxx", ".h++", ".inl", ".ipp")
            documented = any(_text(record[1].find(tag)) for tag in ("briefdescription", "detaileddescription", "inbodydescription"))
            return (header, bool(location.get("declfile")), documented, record[1].get("id"))

        _, canonical = max(records, key=preference)
        canonical_id = canonical.get("id")
        for _, member in records:
            if member.get("id") != canonical_id:
                aliases[member.get("id")] = canonical_id
        for path, member in records:
            if member.get("id") == canonical_id:
                for _, source in records:
                    _merge_descriptions(member, source)
                changed.add(path)

    def remap(parent):
        modified = False
        for child in list(parent):
            if ((child.tag == "memberdef" and child.get("id") in aliases)
                    or (child.tag == "member" and child.get("refid") in aliases)):
                parent.remove(child)
                modified = True
                continue
            for attribute in ("id", "refid"):
                identifier = child.get(attribute)
                if identifier in aliases:
                    child.set(attribute, aliases[identifier])
                    modified = True
            modified = remap(child) or modified
        return modified

    if aliases:
        for path, tree in trees.items():
            if path in excluded_files:
                continue
            if remap(tree.getroot()):
                changed.add(path)
    return changed


def _deduplicate_aliases(trees):
    """Doxygen may emit an unqualified alias for the same qualified class."""
    classes = []
    location_fields = ("file", "line", "column", "bodyfile", "bodystart", "bodyend")
    for path, tree in trees.items():
        compound = tree.getroot().find("compounddef")
        if compound is None or compound.get("kind") not in ("class", "struct", "union"):
            continue
        location = compound.find("location")
        if location is None or not all(location.get(field) for field in location_fields):
            continue
        refid = compound.get("id")
        members = tuple(sorted(
            (member.get("kind"), member.get("id", "").removeprefix(refid),
             member.findtext("definition", "").strip(), member.findtext("argsstring", ""))
            for member in compound.findall(".//memberdef")))
        # Location alone cannot distinguish a template from its instantiations.
        if not members or any(not member[2] for member in members):
            continue
        signature = (compound.get("kind"), tuple(location.get(field) for field in location_fields), members)
        classes.append((path, refid, compound.findtext("compoundname", ""), signature))
    aliases = {}
    removed_files = set()
    for path, refid, name, signature in classes:
        if "::" in name:
            continue
        qualified = [(other_id, other_name) for _, other_id, other_name, other_signature in classes
                     if "::" in other_name and other_name.rsplit("::", 1)[-1] == name
                     and other_signature == signature]
        if qualified:
            canonical, _ = max(qualified, key=lambda item: (item[1].count("::"), item[1]))
            aliases[refid] = canonical
            removed_files.add(path)

    def canonical_id(identifier):
        for alias, canonical in aliases.items():
            if identifier == alias or identifier.startswith(alias + "_1"):
                return canonical + identifier[len(alias):]
        return identifier

    changed_files = set()
    for path, tree in trees.items():
        if path in removed_files:
            continue
        root = tree.getroot()
        # Remove the alias entry before remapping, otherwise index.xml has two copies.
        for compound in list(root.findall("compound")):
            if compound.get("refid") in aliases:
                root.remove(compound)
                changed_files.add(path)
        for parent in root.iter():
            for child in list(parent):
                if child.tag == "innerclass" and child.get("refid") in aliases:
                    parent.remove(child)
                    changed_files.add(path)
        for node in root.iter():
            for attribute in ("id", "refid"):
                identifier = node.get(attribute)
                if identifier and canonical_id(identifier) != identifier:
                    node.set(attribute, canonical_id(identifier))
                    changed_files.add(path)
    return removed_files, changed_files


def prune_xml(xml_dir, project_root=None):
    """Filter private entries and duplicate aliases before Exhale reads generated XML."""
    xml_dir = Path(xml_dir)
    trees = {path: ET.parse(path) for path in sorted(xml_dir.glob("*.xml"))}
    excluded_files, changed_files = _deduplicate_aliases(trees)
    if project_root is not None:
        changed_files.update(_normalize_paths(trees, project_root))
    changed_files.update(_deduplicate_functions(trees, excluded_files))
    excluded = set()
    for path, tree in trees.items():
        if path in excluded_files:
            continue
        compound = tree.getroot().find("compounddef")
        if compound is None:
            continue
        location = compound.find("location")
        if (compound.get("prot") == "private"
                or _private_name(compound.findtext("compoundname", ""))
                or location is not None and location.get("file") == "[STL]"):
            excluded_files.add(path)
            excluded.add(compound.get("id"))
            excluded.update(node.get("id") for node in compound.iter()
                            if node.tag in ("memberdef", "enumvalue"))
        for member in compound.findall(".//memberdef"):
            if (member.get("prot") == "private"
                    or _private_name(member.findtext("name", ""))
                    or compound.get("kind") in ("file", "namespace")
                    and member.get("static") == "yes"):
                excluded.add(member.get("id"))
                excluded.update(value.get("id") for value in member.findall("enumvalue"))
            else:
                excluded.update(value.get("id") for value in member.findall("enumvalue")
                                if _private_name(value.findtext("name", "")))
    excluded.discard(None)

    def remove_references(parent):
        changed = False
        for child in list(parent):
            if child.get("id") in excluded or child.get("refid") in excluded:
                if child.tag == "ref":
                    # Keep prose and source listings readable when their link is private.
                    text = "".join(child.itertext()) + (child.tail or "")
                    position = list(parent).index(child)
                    if position:
                        previous = parent[position - 1]
                        previous.tail = (previous.tail or "") + text
                    else:
                        parent.text = (parent.text or "") + text
                parent.remove(child)
                changed = True
            else:
                changed = remove_references(child) or changed
        return changed

    for path, tree in trees.items():
        if path in excluded_files:
            path.unlink()
        elif remove_references(tree.getroot()) or path in changed_files:
            # Exhale's file scanner expects Doxygen's compact self-closing tags.
            xml = ET.tostring(tree.getroot(), encoding="unicode", xml_declaration=True)
            path.write_text(xml.replace(" />", "/>"), encoding="utf-8")
    return excluded


def _raw_anchor(refid):
    # reST labels normalize case and underscores; Doxygen identifiers must stay exact.
    return f'.. raw:: html\n\n   <span id="{escape(refid, quote=True)}"></span>\n\n'


def add_compound_pages(xml_dir, exhale_root):
    """Add header anchors and group descriptions without re-rendering declarations."""
    output = Path(exhale_root.root_directory).resolve()
    groups = []
    aliases = {}
    for path in sorted(Path(xml_dir).glob("*.xml")):
        compound = ET.parse(path).getroot().find("compounddef")
        if compound is None or compound.get("kind") not in ("file", "namespace", "group"):
            continue
        refid = compound.get("id")
        node = exhale_root.node_by_refid.get(refid)
        if node is None:
            continue
        destination = (output / node.file_name).resolve()
        if not destination.is_relative_to(output):
            raise ValueError(f"Exhale page escapes the generated reference: {destination}")
        anchor = _raw_anchor(refid)
        if compound.get("kind") == "group":
            if not any(compound.find(tag) is not None
                       and "".join(compound.find(tag).itertext()).strip()
                       for tag in ("briefdescription", "detaileddescription", "inbodydescription")):
                continue
            title = compound.findtext("title") or compound.findtext("compoundname")
            name = compound.findtext("compoundname")
            destination.write_text(
                anchor + f".. _{node.link_name}:\n\n{title}\n{'=' * len(title)}\n\n"
                + f".. doxygengroup:: {name}\n   :desc-only:\n   :no-link:\n",
                encoding="utf-8",
            )
            node.file_name = destination.relative_to(output).as_posix()
            groups.append(destination.relative_to(output).with_suffix("").as_posix())
        else:
            contents = destination.read_text(encoding="utf-8")
            if anchor not in contents:
                destination.write_text(anchor + contents, encoding="utf-8")
        aliases[refid] = destination
    if groups:
        root_page = Path(exhale_root.full_root_file_path)
        listing = ("\nDocumentation groups\n--------------------\n\n.. toctree::\n"
                   "   :maxdepth: 1\n\n" + "".join(f"   {group}\n" for group in groups))
        contents = root_page.read_text(encoding="utf-8")
        if listing not in contents:
            root_page.write_text(contents + listing, encoding="utf-8")
    return aliases


def _xml_dir(app):
    return Path(app.confdir) / app.config.breathe_projects[app.config.breathe_default_project]


def setup(app):
    app.connect("builder-inited", lambda app: prune_xml(_xml_dir(app), Path(app.confdir).parent), priority=100)
    app.connect("builder-inited", lambda app: add_compound_pages(_xml_dir(app), app.exhale_root), priority=900)

"""Check generated Python/Cython/C++ APIs and local links in the Sphinx site."""

from html.parser import HTMLParser
import json
from pathlib import Path
import sys
from urllib.parse import unquote, urlsplit
from xml.etree import ElementTree as ET


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.anchors = set()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.anchors.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.anchors.add(attrs["name"])
        if tag in ("a", "link") and attrs.get("href"):
            self.links.append(attrs["href"])
        if tag in ("img", "script") and attrs.get("src"):
            self.links.append(attrs["src"])


def broken_links(root):
    root = root.resolve()
    pages = {path: Page(path.read_text(encoding="utf-8")) for path in root.rglob("*.html")}
    errors = []
    for path, page in pages.items():
        for href in page.links:
            url = urlsplit(href)
            if url.scheme or url.netloc:
                continue
            if url.path.startswith("/"):
                target = root / unquote(url.path).lstrip("/")
            elif url.path:
                target = path.parent / unquote(url.path)
            else:
                target = path
            target = target.resolve()
            if target.is_dir():
                target /= "index.html"
            if not target.is_relative_to(root) or not target.is_file():
                errors.append(f"{path.relative_to(root)}: missing file {href}")
            elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].anchors:
                errors.append(f"{path.relative_to(root)}: missing anchor {href}")
    return errors


def coverage_errors(root, entries, label):
    root = root.resolve()
    pages = {}
    errors = []
    for entry in entries:
        path = (root / (entry["docname"] + ".html")).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            errors.append(f"Missing {label} page: {entry['fullname']}")
            continue
        if entry.get("anchor"):
            if path not in pages:
                pages[path] = Page(path.read_text(encoding="utf-8"))
            if entry["anchor"] not in pages[path].anchors:
                errors.append(f"Missing {label} target: {entry['fullname']}")
    return errors


def sphinx_cpp_targets(root, entries):
    root = root.resolve()
    targets = {}
    for path in sorted((root / "generated/cpp").rglob("*.html")):
        for anchor in Page(path.read_text(encoding="utf-8")).anchors:
            targets.setdefault(anchor, []).append(path)
    rendered = []
    errors = []
    for entry in entries:
        refid = entry["refid"]
        matches = targets.get(refid, [])
        if not matches:
            errors.append(f"Missing Sphinx C++ API target: {entry['fullname']} ({refid})")
        elif len(matches) > 1:
            pages = ", ".join(path.relative_to(root).as_posix() for path in matches)
            errors.append(f"Duplicate Sphinx C++ API target: {entry['fullname']} ({refid}): {pages}")
        else:
            rendered.append({
                **entry, "docname": matches[0].relative_to(root).with_suffix("").as_posix(),
                "anchor": refid,
            })
    return rendered, errors


def documented_cpp_symbols(xml_dir):
    entries = {}
    for path in sorted(xml_dir.glob("*.xml")):
        compound = ET.parse(path).getroot().find("compounddef")
        if compound is None or compound.get("prot") == "private":
            continue
        location = compound.find("location")
        if location is not None and location.get("file") == "[STL]":
            continue
        if any(part.startswith("_") for part in compound.findtext("compoundname", "").split("::")):
            continue
        nodes = [compound]
        for member in compound.findall(".//memberdef"):
            if member.get("prot") == "private":
                continue
            if member.findtext("name", "").startswith("_"):
                continue
            if compound.get("kind") in ("file", "namespace") and member.get("static") == "yes":
                continue
            nodes.extend([member, *(value for value in member.findall("enumvalue")
                                    if not value.findtext("name", "").startswith("_"))])
        for node in nodes:
            if not any(node.find(name) is not None and "".join(node.find(name).itertext()).strip()
                       for name in ("briefdescription", "detaileddescription", "inbodydescription")):
                continue
            identifier = node.get("id")
            entries[identifier] = {
                "fullname": node.findtext("definition") or node.findtext("compoundname") or node.findtext("name"),
                "kind": node.get("kind", "enumvalue"),
                "refid": identifier,
            }
    return list(entries.values())


def check(root):
    required = [
        "index.html", "installation.html", "getting_started.html", "generated/cpp/index.html",
        "generated/pele.basinhopping.BasinHopping.html",
        "generated/pele.optimize.Fire.html", "generated/pele.optimize.ModifiedFireCPP.html",
        "generated/pele.optimize.modifiedfire_cpp.html", "generated/pele.optimize.LBFGS_CPP.html",
        "generated/pele.optimize.CVODEBDFOptimizer.run.html",
        "generated/pele.potentials.ATLJCPP.getEnergy.html",
        "generated/pele.potentials.LJ.html", "generated/pele.potentials.HS_WCA.html",
        "generated/pele.potentials.InversePower.html", "generated/pele.potentials.WCA.html",
        "generated/pele.storage.Database.html",
        "generated/pele.potentials._pele.PairwisePotentialInterface.html",
        "generated/source_api/index.html",
    ]
    errors = [f"Missing generated page: {path}" for path in required if not (root / path).is_file()]
    source_manifest = root / "source_api_manifest.json"
    if source_manifest.is_file():
        entries = json.loads(source_manifest.read_text(encoding="utf-8"))
        errors.extend(coverage_errors(root, entries, "source API"))
        print(f"Checked {len(entries)} documented Python/Cython source targets")
    else:
        errors.append("Missing source API coverage manifest")
    xml_dir = Path(__file__).parent / "doxygen/_build/xml"
    cpp_entries = documented_cpp_symbols(xml_dir)
    rendered_cpp = []
    if cpp_entries:
        documented_refids = {entry["refid"] for entry in cpp_entries}
        required_classes = [
            {"fullname": f"pele::{name}", "kind": "class", "refid": f"classpele_1_1{name}"}
            for name in ("BasePotential", "PairwisePotentialInterface", "SimplePairwisePotential", "CellListPotential", "LBFGS", "MODIFIED_FIRE", "CVODEBDFOptimizer", "GenericMixedDescent")
            if f"classpele_1_1{name}" not in documented_refids
        ]
        rendered, target_errors = sphinx_cpp_targets(root, cpp_entries + required_classes)
        errors.extend(target_errors)
        errors.extend(coverage_errors(root, rendered, "Sphinx C++ API"))
        rendered_cpp = [entry for entry in rendered if entry["refid"] in documented_refids]
    else:
        errors.append("Missing documented Doxygen APIs")
    errors.extend(broken_links(root))
    if errors:
        raise SystemExit("\n".join(errors))
    (root / "doxygen_api_manifest.json").write_text(json.dumps(rendered_cpp, indent=2) + "\n", encoding="utf-8")
    print(f"Checked {len(rendered_cpp)} documented native source targets")
    print(f"Generated APIs and local links checked in {root}")


if __name__ == "__main__":
    check(Path(sys.argv[1] if len(sys.argv) > 1 else "_build/html"))

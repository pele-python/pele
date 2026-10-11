"""Check generated APIs and local links in the combined Sphinx/Doxygen site."""

from html.parser import HTMLParser
from pathlib import Path
import sys
from urllib.parse import unquote, urlsplit


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


def check(root):
    required = [
        "index.html", "installation.html", "getting_started.html", "cpp_reference.html",
        "generated/pele.basinhopping.BasinHopping.html",
        "generated/pele.potentials.LJ.html", "generated/pele.storage.Database.html",
        "doxygen/namespacepele.html", "doxygen/annotated.html",
        "doxygen/classpele_1_1CVODEBDFOptimizer.html",
        "doxygen/classpele_1_1GenericMixedDescent.html",
    ]
    errors = [f"Missing generated page: {path}" for path in required if not (root / path).is_file()]
    cpp = root / "cpp_reference.html"
    if cpp.is_file():
        anchors = Page(cpp.read_text(encoding="utf-8")).anchors
        for name in ("BasePotential", "LBFGS", "CVODEBDFOptimizer", "GenericMixedDescent", "get_energy_gradient"):
            if not any(anchor.startswith("_CPPv") and name in anchor for anchor in anchors):
                errors.append(f"Missing Breathe C++ API target: {name}")
    errors.extend(broken_links(root))
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Generated APIs and local links checked in {root}")


if __name__ == "__main__":
    check(Path(sys.argv[1] if len(sys.argv) > 1 else "_build/html"))

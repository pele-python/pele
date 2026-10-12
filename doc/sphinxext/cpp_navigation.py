"""Link topic navigation to the canonical Exhale C++ reference pages."""

from pathlib import Path


def descendants(nodes, root_refid):
    selected = {root_refid}
    if not any(node.refid == root_refid for node in nodes):
        raise ValueError(f"Missing C++ topic root: {root_refid}")
    while True:
        expanded = selected | {
            node.refid for node in nodes
            if any(refid in selected for _, refid, _ in node.base_compounds)
        }
        if expanded == selected:
            return sorted((node for node in nodes if node.refid in selected), key=lambda node: node.name)
        selected = expanded


def generate(app):
    source = Path(app.srcdir).resolve()
    root = app.exhale_root
    output = (source / root.root_directory).resolve()
    root_page = (source / root.full_root_file_path).resolve()
    root_docname = root_page.relative_to(source).with_suffix("").as_posix()
    nodes = [node for node in root.all_nodes if node.kind in ("class", "struct")]
    topic_docnames = {}
    for topic, title, root_refid in (
        ("optimizers", "C++ optimizers", "classpele_1_1GradientOptimizer"),
        ("potentials", "C++ potentials", "classpele_1_1BasePotential"),
    ):
        page = [title, "=" * len(title), ""]
        for node in descendants(nodes, root_refid):
            canonical = (output / node.file_name).resolve()
            if not canonical.is_relative_to(output):
                raise ValueError(f"C++ topic target escapes the generated reference: {canonical}")
            docname = canonical.relative_to(source).with_suffix("").as_posix()
            page.append(f"- :doc:`{node.name} </{docname}>`")
        destination = output / (topic + ".rst")
        destination.write_text("\n".join(page) + "\n", encoding="utf-8")
        topic_docnames[topic] = destination.relative_to(source).with_suffix("").as_posix()
    listing = ("\nBrowse by topic\n---------------\n\n.. toctree::\n"
               "   :maxdepth: 1\n\n   optimizers\n   potentials\n")
    contents = root_page.read_text(encoding="utf-8")
    if listing not in contents:
        root_page.write_text(contents + listing, encoding="utf-8")
    app.config.html_context["cpp_navigation"] = [
        {"label": label, "docname": docname, "anchor": anchor}
        for label, docname, anchor in (
            ("Overview", root_docname, ""),
            ("Optimizers", topic_docnames["optimizers"], ""),
            ("Potentials", topic_docnames["potentials"], ""),
            ("Classes & interfaces", root_docname, "class-hierarchy"),
            ("Namespaces", root_docname, "namespaces"),
            ("Functions", root_docname, "functions"),
            ("Source files", root_docname, "file-hierarchy"),
        )
    ]


def setup(app):
    app.connect("builder-inited", generate, priority=950)

"""Generate documented public APIs from tracked sources without importing pele."""

import ast
from collections import defaultdict
from hashlib import sha256
from io import StringIO
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import textwrap
import tokenize
from urllib.parse import quote
import warnings

from numpydoc.docscrape_sphinx import get_doc_object


def tracked_sources(root):
    result = subprocess.run(
        ["git", "ls-files", "-z", "--", "pele"], cwd=root, check=True, capture_output=True,
    )
    return [root / os.fsdecode(path) for path in result.stdout.split(b"\0") if path]


def declaration(text, line):
    """Keep the actual source signature, including Cython types and GIL clauses."""
    tokens = tokenize.generate_tokens(StringIO(text).readline)
    start = None
    depth = 0
    for token in tokens:
        if start is None:
            if token.start[0] >= line and token.type == tokenize.NAME and token.string in ("def", "class", "cdef", "cpdef", "async"):
                start = token.start
            else:
                continue
        if token.type == tokenize.OP:
            if token.string in "([{":
                depth += 1
            elif token.string in ")]}":
                depth -= 1
            elif token.string == ":" and depth == 0:
                lines = text.splitlines(keepends=True)
                selected = lines[start[0] - 1:token.end[0]]
                selected[-1] = selected[-1][:token.end[1]]
                return textwrap.dedent("".join(selected)).strip()
    raise ValueError(f"No source declaration found at line {line}")


def public_name(name, kind):
    return not name.startswith("_") or kind == "method" and name.startswith("__") and name.endswith("__")


def python_docs(text, source):
    tree = ast.parse(text, filename=str(source))
    records = []
    if ast.get_docstring(tree):
        records.append(("", "module", tree.body[0].lineno, ast.get_docstring(tree), ""))

    def walk(node, scope=(), parent_kind=None):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                kind = "class" if isinstance(child, ast.ClassDef) else "method" if parent_kind == "class" else "function"
                if not public_name(child.name, kind):
                    continue
                name = (*scope, child.name)
                doc = ast.get_docstring(child)
                if doc:
                    records.append((".".join(name), kind, child.lineno, doc, declaration(text, child.lineno)))
                if kind == "class":
                    walk(child, name, kind)
            else:
                walk(child, scope, parent_kind)

    walk(tree)
    return records


def cython_docs(text, source):
    # Parsing stops before Cython's import analysis or native compilation.
    from Cython.Compiler import Errors
    from Cython.Compiler.TreeFragment import StringParseContext, parse_from_strings

    Errors.init_thread()
    tree = parse_from_strings(
        str(source), text, level="module_pxd" if source.suffix == ".pxd" else "module",
        context=StringParseContext(str(source), cpp=True),
    )
    records = []
    if tree.doc:
        records.append(("", "module", 1, str(tree.doc), ""))

    def walk(node, scope=(), parent_kind=None):
        node_type = type(node).__name__
        name = None
        if node_type in ("PyClassDefNode", "CClassDefNode"):
            name = getattr(node, "name", None) or node.class_name
            kind = "class"
        elif node_type == "DefNode":
            name = node.name
            kind = "method" if parent_kind == "class" else "function"
        elif node_type == "CFuncDefNode":
            if not (node.overridable or node.visibility in ("public", "extern") or node.api):
                return
            name = node.declarator.declared_name()
            kind = "method" if parent_kind == "class" else "function"
        if name:
            if not public_name(name, kind):
                return
            scope = (*scope, name)
            parent_kind = kind
            if node.doc:
                records.append((".".join(scope), kind, node.pos[1], str(node.doc), declaration(text, node.pos[1])))
            if kind != "class":
                return
        for attribute in getattr(node, "child_attrs", ()):
            children = getattr(node, attribute, None)
            for child in children if isinstance(children, list) else [children]:
                if child is not None and hasattr(child, "child_attrs"):
                    walk(child, scope, parent_kind)

    walk(tree)
    return records


def module_prose(doc):
    """Replace module API catalogs with the generated source catalog."""
    lines = doc.splitlines()
    result = []
    index = 0
    while index < len(lines):
        line = lines[index]
        directive = re.match(r"^(\s*)\.\. (?:auto\w+|currentmodule|module)::", line)
        if directive:
            indent = len(directive[1])
            index += 1
            while index < len(lines) and (not lines[index].strip() or len(lines[index]) - len(lines[index].lstrip()) > indent):
                index += 1
        else:
            result.append(line)
            index += 1
    return "\n".join(result)


def render_docstring(doc, kind):
    with warnings.catch_warnings(record=True) as issues:
        warnings.simplefilter("always")
        rendered = str(get_doc_object(
            None, what="class" if kind == "class" else "function", doc=doc,
            config={"show_class_members": False, "xref_param_type": False},
        ))
    # shortcut: Numpydoc warnings render verbatim, fix source headings to restore formatting.
    if issues:
        return ".. code-block:: text\n\n" + textwrap.indent(doc, "   ")
    return rendered


def generate(root, output, revision=None):
    root, output = Path(root).resolve(), Path(output)
    revision = revision or os.environ.get("GITHUB_SHA") or "docs-autogeneration"
    manifest = []
    modules = defaultdict(list)
    occurrences = defaultdict(int)
    for source in sorted(tracked_sources(root)):
        relative = source.relative_to(root)
        if source.suffix not in (".py", ".pyx", ".pxd") or "tests" in relative.parts or source.name.startswith("test_") or source.name.endswith("_test.py"):
            continue
        with tokenize.open(source) as stream:
            text = stream.read()
        module = ".".join(relative.with_suffix("").parts).removesuffix(".__init__")
        records = python_docs(text, source) if source.suffix == ".py" else cython_docs(text, source)
        if any(part.startswith("_") for part in module.split(".")) and not any(record[1] != "module" for record in records):
            records = []
        for name, kind, line, doc, definition in records:
            fullname = f"{module}.{name}" if name else module
            key = f"{relative.as_posix()}:{kind}:{fullname}"
            occurrences[key] += 1
            anchor = "source-api-" + sha256(f"{key}:{occurrences[key]}".encode()).hexdigest()[:20]
            entry = {
                "fullname": fullname, "kind": kind, "source": relative.as_posix(), "line": line,
                "docstring": doc, "docname": f"generated/source_api/{module}", "anchor": anchor,
            }
            manifest.append(entry)
            modules[module].append((entry, name, definition, source.suffix))

    output.mkdir(parents=True, exist_ok=True)
    for stale in output.glob("*.rst"):
        if stale.stem not in modules and stale.name != "index.rst":
            stale.unlink()
    index = ["Developer source API", "====================", "", "Public library and developer interfaces are collected from source without importing optional backends. Private helpers and definitions nested inside functions are excluded.", "", ".. toctree::", "   :maxdepth: 1", ""]
    for module, records in sorted(modules.items()):
        index.append(f"   {module}")
        page = [module, "=" * len(module), "", f".. py:currentmodule:: {module}", ""]
        for entry, name, definition, suffix in records:
            page.extend([f".. _{entry['anchor']}:", ""])
            url = f"https://github.com/pele-python/pele/blob/{quote(revision, safe='')}/{quote(entry['source'])}#L{entry['line']}"
            if entry["kind"] == "module":
                prose = module_prose(entry["docstring"])
                page.extend([f"`Module source <{url}>`__", "", prose, ""])
                continue
            doc = render_docstring(entry["docstring"], entry["kind"])
            body = [f"`Source <{url}>`__", "", f".. code-block:: {'python' if suffix == '.py' else 'cython'}", "", textwrap.indent(definition, "   "), "", doc]
            page.extend([f".. py:{entry['kind']}:: {name}", "   :no-index:", "", textwrap.indent("\n".join(body), "   "), ""])
        (output / f"{module}.rst").write_text("\n".join(page), encoding="utf-8")
    (output / "index.rst").write_text("\n".join(index) + "\n", encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def setup(app):
    app.add_config_value("source_api_revision", None, "env")

    def build_sources(app):
        generate(Path(app.srcdir).parent, Path(app.srcdir) / "generated/source_api", app.config.source_api_revision)

    def copy_manifest(app, exception):
        if exception is None and app.statuscode == 0 and app.builder.format == "html":
            shutil.copyfile(Path(app.srcdir) / "generated/source_api/manifest.json", Path(app.outdir) / "source_api_manifest.json")

    app.connect("builder-inited", build_sources, priority=100)
    app.connect("build-finished", copy_manifest)

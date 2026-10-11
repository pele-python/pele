from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_html import broken_links, check


class CheckHtmlTests(unittest.TestCase):
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
            (root / "cpp_reference.html").write_text("Breathe has not rendered any APIs")
            with self.assertRaises(SystemExit) as result:
                check(root)
            self.assertIn("Missing generated page: index.html", str(result.exception))
            self.assertIn("Missing Breathe C++ API target: BasePotential", str(result.exception))


if __name__ == "__main__":
    unittest.main()

import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from check_version import read_versions


class VersionConsistencyTests(unittest.TestCase):
    def test_release_versions_match(self) -> None:
        versions = read_versions()
        self.assertEqual(len(set(versions.values())), 1, versions)


if __name__ == "__main__":
    unittest.main()


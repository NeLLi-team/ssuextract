import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
CLI = REPO / "scripts" / "pipeline_cli.sh"


def read_config(path: Path) -> str:
    command = 'source "$1"; CONFIG_FILE="$2"; read_config_database_path'
    result = subprocess.run(
        ["bash", "-c", command, "bash", str(CLI), str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.rstrip("\n")


def read_profile(*arguments: str) -> str:
    command = 'source "$1"; shift; read_cli_database_profile "$@"'
    result = subprocess.run(
        ["bash", "-c", command, "bash", str(CLI), *arguments],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.rstrip("\n")


class DatabaseConfigTests(unittest.TestCase):
    def test_plain_path_may_contain_database_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "local.config"
            config.write_text("/data/database_path/db\n")
            self.assertEqual(read_config(config), "/data/database_path/db")

    def test_legacy_nextflow_config_is_still_read(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "local.config"
            config.write_text("params {\n    database_path = '/legacy/db'\n}\n")
            self.assertEqual(read_config(config), "/legacy/db")

    def test_database_profile_accepts_separate_value(self) -> None:
        self.assertEqual(read_profile("--database_profile", "img"), "img")

    def test_database_profile_accepts_equals_value(self) -> None:
        self.assertEqual(read_profile("--database_profile=curated"), "curated")


if __name__ == "__main__":
    unittest.main()

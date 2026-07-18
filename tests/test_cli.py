import os
import re
import subprocess
import stat
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


def read_config_profile(path: Path) -> str:
    command = 'source "$1"; CONFIG_FILE="$2"; read_config_database_profile'
    result = subprocess.run(
        ["bash", "-c", command, "bash", str(CLI), str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.rstrip("\n")


def select_profile(selection: str, default: str = "curated") -> subprocess.CompletedProcess:
    profiles = (
        "curated\t1.0.0\t345.5 MiB\tPR2 and SILVA\n"
        "img\t1.0.0\t828.9 MiB\tIMG enhanced"
    )
    command = 'source "$1"; profile_from_selection "$2" "$3" "$4"'
    return subprocess.run(
        ["bash", "-c", command, "bash", str(CLI), selection, default, profiles],
        capture_output=True,
        text=True,
    )


def write_config(path: Path, database_path: str, profile: str) -> subprocess.CompletedProcess:
    command = 'source "$1"; CONFIG_FILE="$2"; write_database_config "$3" "$4"'
    return subprocess.run(
        [
            "bash",
            "-c",
            command,
            "bash",
            str(CLI),
            str(path),
            database_path,
            profile,
        ],
        capture_output=True,
        text=True,
    )


class DatabaseConfigTests(unittest.TestCase):
    def test_pipeline_cli_remains_directly_executable(self) -> None:
        self.assertTrue(CLI.stat().st_mode & stat.S_IXUSR)

    def test_pipeline_cli_defaults_to_python3(self) -> None:
        command = 'unset PYTHON; source "$1"; printf "%s\\n" "$PYTHON"'
        result = subprocess.run(
            ["bash", "-c", command, "bash", str(CLI)],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.stdout.strip(), "python3")

    def test_version_does_not_launch_nextflow(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            nextflow = executable_dir / "nextflow"
            nextflow.write_text(
                "#!/usr/bin/env bash\n"
                "printf 'unexpected nextflow call\\n' >&2\n"
                "exit 99\n"
            )
            nextflow.chmod(0o755)
            result = subprocess.run(
                [
                    "bash",
                    "-c",
                    'source "$1"; run_pipeline --version',
                    "bash",
                    str(CLI),
                ],
                check=True,
                capture_output=True,
                text=True,
                env={
                    **os.environ,
                    "PATH": f"{executable_dir}:{os.environ['PATH']}",
                },
            )
        self.assertEqual(result.stdout.strip(), "1.1.0")
        self.assertEqual(result.stderr, "")

    def test_nextflow_uses_safe_term_when_current_term_is_unsupported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            tput = executable_dir / "tput"
            nextflow = executable_dir / "nextflow"
            tput.write_text("#!/usr/bin/env bash\nexit 1\n")
            nextflow.write_text("#!/usr/bin/env bash\nprintf '%s|%s\\n' \"$TERM\" \"$*\"\n")
            tput.chmod(0o755)
            nextflow.chmod(0o755)
            result = subprocess.run(
                [
                    "bash",
                    "-c",
                    'source "$1"; run_nextflow run workflow.nf --help',
                    "bash",
                    str(CLI),
                ],
                check=True,
                capture_output=True,
                text=True,
                env={
                    **os.environ,
                    "PATH": f"{executable_dir}:{os.environ['PATH']}",
                    "TERM": "xterm-kitty",
                },
            )
        self.assertEqual(result.stdout.strip(), "xterm-256color|run workflow.nf --help")

    def test_nextflow_keeps_supported_term(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            tput = executable_dir / "tput"
            nextflow = executable_dir / "nextflow"
            tput.write_text("#!/usr/bin/env bash\nexit 0\n")
            nextflow.write_text("#!/usr/bin/env bash\nprintf '%s|%s\\n' \"$TERM\" \"$*\"\n")
            tput.chmod(0o755)
            nextflow.chmod(0o755)
            result = subprocess.run(
                [
                    "bash",
                    "-c",
                    'source "$1"; run_nextflow run workflow.nf --help',
                    "bash",
                    str(CLI),
                ],
                check=True,
                capture_output=True,
                text=True,
                env={
                    **os.environ,
                    "PATH": f"{executable_dir}:{os.environ['PATH']}",
                    "TERM": "xterm-kitty",
                },
            )
        self.assertEqual(result.stdout.strip(), "xterm-kitty|run workflow.nf --help")

    def test_pipeline_injects_configured_database_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            tput = executable_dir / "tput"
            nextflow = executable_dir / "nextflow"
            tput.write_text("#!/usr/bin/env bash\nexit 0\n")
            nextflow.write_text("#!/usr/bin/env bash\nprintf '%s\\n' \"$*\"\n")
            tput.chmod(0o755)
            nextflow.chmod(0o755)
            command = (
                'source "$1"; '
                'resolve_database_profile() { printf "curated\\n"; }; '
                'resolve_database_path() { printf "/database\\n"; }; '
                'write_database_config() { :; }; ensure_database() { :; }; '
                'check_database_update() { :; }; '
                'run_pipeline --querydir /queries --outdir /output'
            )
            result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI)],
                check=True,
                capture_output=True,
                text=True,
                env={
                    **os.environ,
                    "PATH": f"{executable_dir}:{os.environ['PATH']}",
                    "TERM": "xterm-kitty",
                },
            )
        self.assertEqual(
            result.stdout.strip(),
            f"run {REPO / 'main.nf'} --database_path /database "
            "--database_profile curated --querydir /queries --outdir /output",
        )

    def test_public_commands_do_not_use_a_pixi_argument_separator(self) -> None:
        public_docs = [
            path
            for path in sorted((REPO / "docs").rglob("*.md"))
            if "handoffs" not in path.parts
        ]
        paths = [
            REPO / "README.md",
            REPO / "main.nf",
            *public_docs,
            *sorted((REPO / "scripts").glob("*.py")),
            *sorted((REPO / "scripts").glob("*.sh")),
        ]
        text = "\n".join(path.read_text() for path in paths)
        self.assertNotIn("pixi run setup -- --", text)
        self.assertNotIn("pixi run ssuextract -- --", text)
        self.assertIsNone(
            re.search(r"^pixi run (?:setup|ssuextract).*\\$", text, re.MULTILINE)
        )

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

    def test_current_config_stores_path_and_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "local.config"
            config.write_text(
                "database_path=/data/ssuextract\ndatabase_profile=img\n"
            )
            self.assertEqual(read_config(config), "/data/ssuextract")
            self.assertEqual(read_config_profile(config), "img")

    def test_legacy_nextflow_profile_is_still_read(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "local.config"
            config.write_text("params {\n    database_profile = 'img'\n}\n")
            self.assertEqual(read_config_profile(config), "img")

    def test_database_profile_accepts_separate_value(self) -> None:
        self.assertEqual(read_profile("--database_profile", "img"), "img")

    def test_database_profile_accepts_equals_value(self) -> None:
        self.assertEqual(read_profile("--database_profile=curated"), "curated")

    def test_profile_selection_accepts_number_name_and_default(self) -> None:
        self.assertEqual(select_profile("2").stdout.strip(), "img")
        self.assertEqual(select_profile("img").stdout.strip(), "img")
        self.assertEqual(select_profile("").stdout.strip(), "curated")

    def test_profile_selection_rejects_unknown_value(self) -> None:
        result = select_profile("unknown")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_profile_selection_rejects_out_of_range_number(self) -> None:
        result = select_profile("99")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_config_writer_rejects_unsafe_values(self) -> None:
        cases = (("", "curated"), ("/path\nother=value", "curated"), ("/path", "bad/name"))
        with tempfile.TemporaryDirectory() as tmp:
            for number, (database_path, profile) in enumerate(cases):
                with self.subTest(database_path=database_path, profile=profile):
                    config = Path(tmp) / f"local-{number}.config"
                    result = write_config(config, database_path, profile)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertFalse(config.exists())

    def test_profile_prompt_fails_if_profiles_cannot_be_listed(self) -> None:
        command = 'source "$1"; PYTHON=false; prompt_database_profile curated'
        result = subprocess.run(
            ["bash", "-c", command, "bash", str(CLI)],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Could not list database profiles", result.stderr)

    def test_profile_prompt_labels_saved_default_and_keeps_both_choices(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fake_python = Path(tmp) / "fake-python"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "printf 'curated\\t1.0.1\\t345.4 MiB\\tPR2 and SILVA\\n'\n"
                "printf 'img\\t1.0.1\\t828.8 MiB\\tIMG enhanced\\n'\n"
            )
            fake_python.chmod(0o755)
            command = 'source "$1"; PYTHON="$2"; prompt_database_profile img'
            result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(fake_python)],
                input="1\n",
                check=True,
                capture_output=True,
                text=True,
            )
            default_result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(fake_python)],
                input="\n",
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.stdout, "curated\n")
        self.assertEqual(default_result.stdout, "img\n")
        self.assertIn("  1) curated v1.0.1", result.stderr)
        self.assertIn("  2) img v1.0.1", result.stderr)
        self.assertIn(
            "Database profile (1-2) [default: 2 (img)]: ", result.stderr
        )

    def test_run_update_notice_is_written_to_stderr(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fake_python = Path(tmp) / "fake-python"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$2\" == version ]]; then\n"
                "  echo 1.0.0\n"
                "elif [[ \"$2\" == check-update ]]; then\n"
                "  echo \"Database profile 'curated': installed v1.0.0; latest v1.0.0.\"\n"
                "fi\n"
            )
            fake_python.chmod(0o755)
            command = (
                'source "$1"; PYTHON="$2"; '
                'check_database_update /database curated'
            )
            result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(fake_python)],
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.stdout, "")
        self.assertIn("Checking Zenodo for database updates", result.stderr)
        self.assertIn("installed v1.0.0; latest v1.0.0", result.stderr)


if __name__ == "__main__":
    unittest.main()

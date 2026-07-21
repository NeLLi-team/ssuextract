import os
import re
import subprocess
import stat
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
CLI = REPO / "scripts" / "pipeline_cli.sh"
DATABASE_MANAGER = REPO / "scripts" / "database_manager.py"


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
        self.assertEqual(result.stdout.strip(), "1.2.1")
        self.assertEqual(result.stderr, "")

    def test_example_uses_both_bundled_assemblies(self) -> None:
        command = 'source "$1"; printf "%s\\n" "${SMOKE_FASTAS[@]}"'
        result = subprocess.run(
            ["bash", "-c", command, "bash", str(CLI)],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            result.stdout.splitlines(),
            [
                str(REPO / "data" / "example" / "LKH462_P08_Rh.fna"),
                str(REPO / "data" / "example" / "LKH565_P11_Ci.fna"),
            ],
        )

    def test_example_version_does_not_require_a_database(self) -> None:
        command = 'source "$1"; run_smoke --version'
        result = subprocess.run(
            ["bash", "-c", command, "bash", str(CLI)],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.stdout.strip(), "1.2.1")
        self.assertEqual(result.stderr, "")

    def test_example_default_wiring_persists_curated_and_uses_two_threads(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            temporary = Path(tmp)
            first = temporary / "first.fna"
            second = temporary / "second.fna"
            fake_python = temporary / "fake-python"
            captured_config = temporary / "config.txt"
            captured_arguments = temporary / "arguments.txt"
            first.write_text(">first\nACGT\n")
            second.write_text(">second\nTGCA\n")
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$1\" == *database_manager.py && \"$2\" == version ]]; then\n"
                "  printf '1.0.1\\n'\n"
                "elif [[ \"$1\" == *validate_example_output.py ]]; then\n"
                "  exit 0\n"
                "else\n"
                "  exit 99\n"
                "fi\n"
            )
            fake_python.chmod(0o755)
            command = (
                'source "$1"; '
                'SMOKE_FASTAS=("$2" "$3"); PYTHON="$4"; '
                'CAPTURED_CONFIG="$5"; CAPTURED_ARGUMENTS="$6"; '
                'resolve_database_profile() { printf "curated\\n"; }; '
                'resolve_database_path() { printf "/database\\n"; }; '
                'write_database_config() { printf "%s|%s\\n" "$1" "$2" '
                '> "$CAPTURED_CONFIG"; }; '
                'run_pipeline() { printf "%s\\n" "$@" '
                '> "$CAPTURED_ARGUMENTS"; }; '
                'run_smoke'
            )
            subprocess.run(
                [
                    "bash",
                    "-c",
                    command,
                    "bash",
                    str(CLI),
                    str(first),
                    str(second),
                    str(fake_python),
                    str(captured_config),
                    str(captured_arguments),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            arguments = captured_arguments.read_text().splitlines()
            config = captured_config.read_text().strip()
        self.assertEqual(config, "/database|curated")
        self.assertEqual(arguments[arguments.index("--threads_per_job") + 1], "2")
        self.assertNotIn("--database_path", arguments)
        self.assertNotIn("--database_profile", arguments)

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

    def test_pipeline_disables_update_prompt_for_explicit_database_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            tput = executable_dir / "tput"
            nextflow = executable_dir / "nextflow"
            captured = executable_dir / "update-prompt.txt"
            tput.write_text("#!/usr/bin/env bash\nexit 0\n")
            nextflow.write_text("#!/usr/bin/env bash\nexit 0\n")
            tput.chmod(0o755)
            nextflow.chmod(0o755)
            command = (
                'source "$1"; CAPTURED="$2"; '
                'ensure_database() { :; }; '
                'check_database_update() { printf "%s|%s\\n" "$3" "$5" '
                '> "$CAPTURED"; }; '
                'run_pipeline --database_path /custom --database_profile img '
                '--querydir /queries --outdir /output'
            )
            subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(captured)],
                check=True,
                capture_output=True,
                text=True,
                env={
                    **os.environ,
                    "PATH": f"{executable_dir}:{os.environ['PATH']}",
                    "TERM": "xterm-kitty",
                },
            )
            captured_value = captured.read_text().strip()
        self.assertEqual(captured_value, "0|1")

    def test_pipeline_allows_update_prompt_for_saved_database_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            tput = executable_dir / "tput"
            nextflow = executable_dir / "nextflow"
            captured = executable_dir / "update-prompt.txt"
            tput.write_text("#!/usr/bin/env bash\nexit 0\n")
            nextflow.write_text("#!/usr/bin/env bash\nexit 0\n")
            tput.chmod(0o755)
            nextflow.chmod(0o755)
            command = (
                'source "$1"; CAPTURED="$2"; '
                'resolve_database_profile() { printf "img\\n"; }; '
                'resolve_database_path() { printf "/managed\\n"; }; '
                'write_database_config() { :; }; ensure_database() { :; }; '
                'check_database_update() { printf "%s|%s\\n" "$3" "$5" '
                '> "$CAPTURED"; }; '
                'run_pipeline --querydir /queries --outdir /output'
            )
            subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(captured)],
                check=True,
                capture_output=True,
                text=True,
                env={
                    **os.environ,
                    "PATH": f"{executable_dir}:{os.environ['PATH']}",
                    "TERM": "xterm-kitty",
                },
            )
            captured_value = captured.read_text().strip()
        self.assertEqual(captured_value, "1|0")

    def test_pipeline_continues_after_noninteractive_database_notices(self) -> None:
        cases = (
            (
                "unavailable\t1.0.1\t-\toffline",
                "Database update check unavailable: offline",
            ),
            (
                "update_available\t1.0.1\t1.0.2\t-",
                "Run: pixi run setup --database_profile curated --update",
            ),
        )
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            tput = executable_dir / "tput"
            nextflow = executable_dir / "nextflow"
            fake_python = executable_dir / "fake-python"
            tput.write_text("#!/usr/bin/env bash\nexit 0\n")
            nextflow.write_text("#!/usr/bin/env bash\nprintf '%s\\n' \"$*\"\n")
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$2\" == version ]]; then\n"
                "  echo 1.0.1\n"
                "elif [[ \"$2\" == check-update ]]; then\n"
                "  printf '%s\\n' \"$UPDATE_STATUS\"\n"
                "fi\n"
            )
            for executable in (tput, nextflow, fake_python):
                executable.chmod(0o755)
            command = (
                'source "$1"; PYTHON="$2"; '
                'resolve_database_profile() { printf "curated\\n"; }; '
                'resolve_database_path() { printf "/database\\n"; }; '
                'write_database_config() { :; }; ensure_database() { :; }; '
                'run_pipeline --querydir /queries --outdir /output'
            )
            for update_status, expected_notice in cases:
                with self.subTest(update_status=update_status):
                    result = subprocess.run(
                        [
                            "bash",
                            "-c",
                            command,
                            "bash",
                            str(CLI),
                            str(fake_python),
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                        env={
                            **os.environ,
                            "PATH": f"{executable_dir}:{os.environ['PATH']}",
                            "TERM": "xterm-kitty",
                            "UPDATE_STATUS": update_status,
                        },
                    )
                    self.assertIn(f"run {REPO / 'main.nf'}", result.stdout)
                    self.assertIn(expected_notice, result.stderr)
                    self.assertNotIn("Update database now?", result.stderr)

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
                "  printf 'current\\t1.0.0\\t1.0.0\\t-\\n'\n"
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
        self.assertIn("Database profile curated is current at v1.0.0", result.stderr)

    def test_run_update_prompt_installs_before_pipeline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            temporary = Path(tmp)
            fake_python = temporary / "fake-python"
            captured = temporary / "install-arguments.txt"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "case \"$2\" in\n"
                "  version) echo 1.0.2 ;;\n"
                "  check-update) printf 'update_available\\t1.0.1\\t1.0.2\\t-\\n' ;;\n"
                "  install) printf '%s\\n' \"$*\" > \"$CAPTURED_INSTALL\" ;;\n"
                "  validate) exit 0 ;;\n"
                "  *) exit 99 ;;\n"
                "esac\n"
            )
            fake_python.chmod(0o755)
            command = (
                'source "$1"; PYTHON="$2"; '
                'can_prompt_for_database_update() { return 0; }; '
                'check_database_update /database img 1'
            )
            result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(fake_python)],
                input="y\n",
                check=True,
                capture_output=True,
                text=True,
                env={**os.environ, "CAPTURED_INSTALL": str(captured)},
            )
            captured_arguments = captured.read_text().strip()
        self.assertEqual(
            captured_arguments,
            f"{DATABASE_MANAGER} install --root /database --profile img --latest "
            "--require-version 1.0.2 --force",
        )
        self.assertIn(
            "Database update available for img: v1.0.1 -> v1.0.2",
            result.stderr,
        )
        self.assertIn("Update database now? [y/N]:", result.stderr)
        self.assertIn("Database profile img v1.0.2 ready at /database", result.stderr)

    def test_run_update_prompt_decline_keeps_installed_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            temporary = Path(tmp)
            fake_python = temporary / "fake-python"
            captured = temporary / "install-arguments.txt"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$2\" == version ]]; then\n"
                "  echo 1.0.1\n"
                "elif [[ \"$2\" == check-update ]]; then\n"
                "  printf 'update_available\\t1.0.1\\t1.0.2\\t-\\n'\n"
                "elif [[ \"$2\" == install ]]; then\n"
                "  printf '%s\\n' \"$*\" > \"$CAPTURED_INSTALL\"\n"
                "fi\n"
            )
            fake_python.chmod(0o755)
            command = (
                'source "$1"; PYTHON="$2"; '
                'can_prompt_for_database_update() { return 0; }; '
                'check_database_update /database curated 1'
            )
            result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(fake_python)],
                input="n\n",
                check=True,
                capture_output=True,
                text=True,
                env={**os.environ, "CAPTURED_INSTALL": str(captured)},
            )
            install_was_called = captured.exists()
        self.assertFalse(install_was_called)
        self.assertIn("Update database now? [y/N]:", result.stderr)
        self.assertIn(
            "pixi run setup --database_profile curated --update",
            result.stderr,
        )

    def test_run_update_notice_includes_explicit_database_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fake_python = Path(tmp) / "fake-python"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$2\" == version ]]; then\n"
                "  echo 1.0.1\n"
                "elif [[ \"$2\" == check-update ]]; then\n"
                "  printf 'update_available\\t1.0.1\\t1.0.2\\t-\\n'\n"
                "fi\n"
            )
            fake_python.chmod(0o755)
            command = (
                'source "$1"; PYTHON="$2"; '
                'check_database_update "/custom database" curated 0 0 1'
            )
            result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(fake_python)],
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertIn(
            "--database_path /custom\\ database --database_profile curated --update",
            result.stderr,
        )

    def test_failed_run_update_stops_before_pipeline_with_old_profile_intact(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fake_python = Path(tmp) / "fake-python"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "case \"$2\" in\n"
                "  version) echo 1.0.1 ;;\n"
                "  check-update) printf 'update_available\\t1.0.1\\t1.0.2\\t-\\n' ;;\n"
                "  install) exit 2 ;;\n"
                "  *) exit 99 ;;\n"
                "esac\n"
            )
            fake_python.chmod(0o755)
            command = (
                'source "$1"; PYTHON="$2"; '
                'can_prompt_for_database_update() { return 0; }; '
                'check_database_update /database img 1'
            )
            result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(fake_python)],
                input="y\n",
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "Database update failed; profile img remains installed and the pipeline was not started.",
            result.stderr,
        )

    def test_setup_update_forces_shared_update_install_path(self) -> None:
        command = (
            'source "$1"; '
            'resolve_setup_database_profile() { printf "img\\n"; }; '
            'resolve_database_path() { printf "/database\\n"; }; '
            'write_database_config() { :; }; database_ready() { return 0; }; '
            'ensure_database() { :; }; database_version() { printf "1.0.1\\n"; }; '
            'check_database_update() { printf "%s|%s|%s|%s\\n" "$@"; }; '
            'setup_database --database_profile img --update'
        )
        result = subprocess.run(
            ["bash", "-c", command, "bash", str(CLI)],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.stdout.strip(), "/database|img|1|1")


if __name__ == "__main__":
    unittest.main()

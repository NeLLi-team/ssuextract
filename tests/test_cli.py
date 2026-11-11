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


def read_config_failure(path: Path, reader: str) -> subprocess.CompletedProcess:
    command = 'source "$1"; CONFIG_FILE="$2"; "$3"'
    return subprocess.run(
        ["bash", "-c", command, "bash", str(CLI), str(path), reader],
        check=False,
        capture_output=True,
        text=True,
    )


FAKE_MANAGER = """#!/usr/bin/env bash
printf '%s\\n' "$*" >> "$CALLS"
if [[ "$1" == *validate_example_output.py ]]; then exit 0; fi
case "$2" in
  validate)
    if [[ -e "$INSTALLED" ]]; then exit 0; fi
    echo "database-manager: $VALIDATE_MESSAGE" >&2
    exit "$VALIDATE_STATUS" ;;
  install) touch "$INSTALLED" ;;
  version) echo 1.0.1 ;;
  profiles)
    printf 'curated\\t1.0.1\\t345.4 MiB\\tPR2 and SILVA\\n'
    printf 'img\\t1.0.1\\t828.8 MiB\\tIMG enhanced\\n' ;;
  check-update) printf '%s\\n' "$UPDATE_STATUS" ;;
  *) exit 99 ;;
esac
"""
VALIDATE_MESSAGE = "validation exploded"


class Sandbox:
    """A wrapper run with a fake database manager and a fake nextflow."""

    def __init__(
        self,
        root: Path,
        validate_status: int = 0,
        installed: bool = False,
        update_status: str = "current\t1.0.1\t1.0.1\t-",
    ) -> None:
        bin_dir = root / "bin"
        bin_dir.mkdir()
        self.python = root / "fake-python"
        self.python.write_text(FAKE_MANAGER)
        self.python.chmod(0o755)
        (bin_dir / "nextflow").write_text(
            '#!/usr/bin/env bash\nprintf \'NEXTFLOW %s\\n\' "$*" >> "$CALLS"\n'
        )
        (bin_dir / "tput").write_text("#!/usr/bin/env bash\nexit 0\n")
        for executable in bin_dir.iterdir():
            executable.chmod(0o755)
        self.calls = root / "calls.txt"
        self.config_dir = root / "config"
        self.config_dir.mkdir()
        self.config = self.config_dir / "local.config"
        self.env = {
            **os.environ,
            "PATH": f"{bin_dir}:{os.environ['PATH']}",
            "TERM": "xterm-kitty",
            "CALLS": str(self.calls),
            "INSTALLED": str(root / "installed"),
            "VALIDATE_STATUS": str(validate_status),
            "VALIDATE_MESSAGE": VALIDATE_MESSAGE,
            "UPDATE_STATUS": update_status,
        }
        if installed:
            (root / "installed").touch()

    def run(self, command: str, *arguments: str) -> subprocess.CompletedProcess:
        script = (
            'source "$1"; PYTHON="$2"; CONFIG_FILE="$3"; DEFAULT_DB_DIR=/database; '
            "shift 3; " + command
        )
        return subprocess.run(
            [
                "bash",
                "-c",
                script,
                "bash",
                str(CLI),
                str(self.python),
                str(self.config),
                *arguments,
            ],
            check=False,
            capture_output=True,
            text=True,
            stdin=subprocess.DEVNULL,
            env=self.env,
        )

    def manager_calls(self) -> list[str]:
        if not self.calls.exists():
            return []
        return self.calls.read_text().splitlines()


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
    def test_array_expansions_are_guarded_for_bash32_nounset(self) -> None:
        source = CLI.read_text(encoding="utf-8")
        guarded = r'\$\{([A-Za-z_]\w*)\[@\]\+"\$\{\1\[@\]\}"\}'
        unguarded = re.findall(r'"\$\{\w+\[@\]\}"', re.sub(guarded, "", source))
        self.assertEqual(unguarded, [], "unsafe empty-array expansions under bash 3.2")

    def test_setup_handles_empty_validation_and_update_arrays_under_nounset(
        self,
    ) -> None:
        for installed in (False, True):
            with (
                self.subTest(installed=installed),
                tempfile.TemporaryDirectory() as tmp,
            ):
                sandbox = Sandbox(Path(tmp), validate_status=3, installed=installed)
                completed = sandbox.run("set -u; setup_database --database_profile img")
                self.assertEqual(completed.returncode, 0, completed.stderr)
                calls = sandbox.manager_calls()
                self.assertTrue(any(" validate " in call for call in calls))
                self.assertTrue(any(" check-update " in call for call in calls))

    def test_workflow_single_paths_disable_glob_expansion(self) -> None:
        source = (REPO / "main.nf").read_text(encoding="utf-8")
        calls = re.findall(r"\bfile\(([^\n]*)\)", source)
        self.assertTrue(calls)
        for arguments in calls:
            with self.subTest(arguments=arguments):
                self.assertRegex(arguments, r"\bglob:\s*false\b")

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
        self.assertEqual(result.stdout.strip(), "1.2.2")
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
        self.assertEqual(result.stdout.strip(), "1.2.2")
        self.assertEqual(result.stderr, "")

    def test_example_default_wiring_uses_two_threads_and_saved_selection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            temporary = Path(tmp)
            first = temporary / "first.fna"
            second = temporary / "second.fna"
            fake_python = temporary / "fake-python"
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
                'CAPTURED_ARGUMENTS="$5"; '
                'resolve_database_profile() { printf "curated\\n"; }; '
                'resolve_database_path() { printf "/database\\n"; }; '
                "run_pipeline() { RUN_DATABASE_PATH=/database; "
                "RUN_DATABASE_PROFILE=curated; "
                'printf "%s\\n" "$@" '
                '> "$CAPTURED_ARGUMENTS"; }; '
                "run_smoke"
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
                    str(captured_arguments),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            arguments = captured_arguments.read_text().splitlines()
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
                "write_database_config() { :; }; ensure_database() { :; }; "
                "require_known_profile() { :; }; "
                "check_database_update() { :; }; "
                "run_pipeline --query /queries --outdir /output"
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
            "--database_profile curated --query /queries --outdir /output",
        )

    def test_query_short_alias_preserves_paths_and_reaches_nextflow_as_query(self) -> None:
        cases = (
            (["-q", "/queries with spaces"], "/queries with spaces"),
            (["-q=/queries-equals"], "/queries-equals"),
        )
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            tput = executable_dir / "tput"
            nextflow = executable_dir / "nextflow"
            tput.write_text("#!/usr/bin/env bash\nexit 0\n")
            nextflow.write_text("#!/usr/bin/env bash\nprintf '<%s>\\n' \"$@\"\n")
            tput.chmod(0o755)
            nextflow.chmod(0o755)
            command = (
                'source "$1"; shift; '
                'resolve_database_profile() { printf "curated\\n"; }; '
                'resolve_database_path() { printf "/database\\n"; }; '
                "write_database_config() { :; }; ensure_database() { :; }; "
                "require_known_profile() { :; }; "
                'check_database_update() { :; }; run_pipeline "$@"'
            )
            for arguments, query_path in cases:
                with self.subTest(arguments=arguments):
                    result = subprocess.run(
                        ["bash", "-c", command, "bash", str(CLI), *arguments],
                        check=True,
                        capture_output=True,
                        text=True,
                        env={
                            **os.environ,
                            "PATH": f"{executable_dir}:{os.environ['PATH']}",
                            "TERM": "xterm-kitty",
                        },
                    )
                    forwarded = result.stdout.splitlines()
                    self.assertIn("<--query>", forwarded)
                    self.assertIn(f"<{query_path}>", forwarded)
                    self.assertNotIn("<-q>", forwarded)
                    self.assertFalse(
                        any(line.startswith("<-q=") for line in forwarded)
                    )

    def test_query_short_alias_requires_a_value_before_database_setup(self) -> None:
        for arguments in (
            ("-q",),
            ("-q", ""),
            ("-q=",),
            ("-q", "--outdir", "/output"),
        ):
            with self.subTest(arguments=arguments):
                command = (
                    'source "$1"; shift; '
                    'resolve_database_profile() { '
                    'printf "database setup reached\\n" >&2; return 99; }; '
                    'run_pipeline "$@"'
                )
                result = subprocess.run(
                    ["bash", "-c", command, "bash", str(CLI), *arguments],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 2)
                self.assertEqual(
                    result.stderr,
                    "-q requires a FASTA file or directory path.\n",
                )

    def test_no_argument_run_preserves_default_nextflow_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            executable_dir = Path(tmp)
            tput = executable_dir / "tput"
            nextflow = executable_dir / "nextflow"
            tput.write_text("#!/usr/bin/env bash\nexit 0\n")
            nextflow.write_text("#!/usr/bin/env bash\nprintf '<%s>\\n' \"$@\"\n")
            tput.chmod(0o755)
            nextflow.chmod(0o755)
            command = (
                'source "$1"; '
                'resolve_database_profile() { printf "curated\\n"; }; '
                'resolve_database_path() { printf "/database\\n"; }; '
                "write_database_config() { :; }; ensure_database() { :; }; "
                "require_known_profile() { :; }; "
                "check_database_update() { :; }; run_pipeline"
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
            result.stdout.splitlines(),
            [
                "<run>",
                f"<{REPO / 'main.nf'}>",
                "<--database_path>",
                "</database>",
                "<--database_profile>",
                "<curated>",
            ],
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
                "ensure_database() { :; }; require_known_profile() { :; }; "
                'check_database_update() { printf "%s|%s\\n" "$3" "$5" '
                '> "$CAPTURED"; }; '
                "run_pipeline --database_path /custom --database_profile img "
                "--query /queries --outdir /output"
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
                "write_database_config() { :; }; ensure_database() { :; }; "
                "require_known_profile() { :; }; "
                'check_database_update() { printf "%s|%s\\n" "$3" "$5" '
                '> "$CAPTURED"; }; '
                "run_pipeline --query /queries --outdir /output"
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
                "write_database_config() { :; }; ensure_database() { :; }; "
                "require_known_profile() { :; }; "
                "run_pipeline --query /queries --outdir /output"
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
                        stdin=subprocess.DEVNULL,
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

    def test_config_bare_first_line_is_not_read_as_a_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "local.config"
            config.write_text("/data/database_path/db\n")
            result = read_config_failure(config, "read_config_database_path")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_config_nextflow_params_block_is_not_read(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "local.config"
            config.write_text(
                "params {\n    database_path = '/old/db'\n"
                "    database_profile = 'img'\n}\n"
            )
            path_result = read_config_failure(config, "read_config_database_path")
            profile_result = read_config_failure(config, "read_config_database_profile")
        self.assertNotEqual(path_result.returncode, 0)
        self.assertEqual(path_result.stdout, "")
        self.assertNotEqual(profile_result.returncode, 0)
        self.assertEqual(profile_result.stdout, "")

    def test_config_with_only_a_profile_yields_no_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "local.config"
            config.write_text("database_profile=img\n")
            result = read_config_failure(config, "read_config_database_path")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_current_config_stores_path_and_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "local.config"
            config.write_text(
                "database_path=/data/ssuextract\ndatabase_profile=img\n"
            )
            self.assertEqual(read_config(config), "/data/ssuextract")
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

    def test_cached_update_does_not_claim_a_zenodo_request(self) -> None:
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
                "check_database_update /database curated 0 0 0 --cached"
            )
            result = subprocess.run(
                ["bash", "-c", command, "bash", str(CLI), str(fake_python)],
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.stdout, "")
        self.assertNotIn("Checking Zenodo for database updates", result.stderr)
        self.assertIn("Database profile curated is current at v1.0.0", result.stderr)

    def test_run_update_prompt_installs_before_pipeline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(
                Path(tmp),
                installed=True,
                update_status="update_available\t1.0.1\t1.0.2\t-",
            )
            result = sandbox.run(
                "can_prompt_for_database_update() { return 0; }; "
                'read() { if [[ "$*" == "-r answer" ]]; then answer=y; '
                'else builtin read "$@"; fi; }; '
                "run_pipeline --database_profile img --query /q"
            )
            calls = sandbox.manager_calls()
        self.assertEqual(result.returncode, 0, result.stderr)
        install_index = next(i for i, call in enumerate(calls) if " install " in call)
        nextflow_index = next(
            i for i, call in enumerate(calls) if call.startswith("NEXTFLOW ")
        )
        self.assertLess(install_index, nextflow_index)
        self.assertIn("--require-version 1.0.2 --force", calls[install_index])

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
            root = Path(tmp)
            sandbox = Sandbox(
                root, installed=True, update_status="update_available\t1.0.1\t1.0.2\t-"
            )
            installed = root / "installed"
            installed.write_text("old profile")
            sandbox.python.write_text(
                FAKE_MANAGER.replace('install) touch "$INSTALLED"', "install) exit 2")
            )
            result = sandbox.run(
                "can_prompt_for_database_update() { return 0; }; "
                'read() { if [[ "$*" == "-r answer" ]]; then answer=y; '
                'else builtin read "$@"; fi; }; '
                "run_pipeline --database_profile img --query /q || exit $?"
            )
            self.assertEqual(installed.read_text(), "old profile")
            calls = sandbox.manager_calls()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(any(call.startswith("NEXTFLOW ") for call in calls))
        self.assertTrue(any(" install " in call for call in calls))
        self.assertIn("Database update failed", result.stderr)

    def test_setup_update_forces_shared_update_install_path(self) -> None:
        command = (
            'source "$1"; '
            'resolve_setup_database_profile() { printf "img\\n"; }; '
            'resolve_database_path() { printf "/database\\n"; }; '
            "write_database_config() { :; }; require_known_profile() { :; }; "
            "ensure_database() { :; }; "
            'check_database_update() { printf "%s|%s|%s|%s\\n" "$@"; }; '
            "setup_database --database_profile img --update"
        )
        result = subprocess.run(
            ["bash", "-c", command, "bash", str(CLI)],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.stdout.strip(), "/database|img|1|1")

    def test_profile_selection_rejects_zero_without_error_output(self) -> None:
        result = select_profile("0")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertEqual(result.stderr, "")

    def test_path_prompt_expands_home_and_stores_absolute_paths(self) -> None:
        cases = (
            ("~/db\n", "/home/fake/db"),
            ("~\n", "/home/fake"),
            ("\n", "/default"),
        )
        with tempfile.TemporaryDirectory() as tmp:
            for answer, expected in cases:
                with self.subTest(answer=answer):
                    result = subprocess.run(
                        [
                            "bash",
                            "-c",
                            'source "$1"; prompt_database_path /default',
                            "bash",
                            str(CLI),
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                        input=answer,
                        env={**os.environ, "HOME": "/home/fake"},
                    )
                    self.assertEqual(result.stdout, expected + "\n")
            relative = subprocess.run(
                [
                    "bash",
                    "-c",
                    'source "$1"; prompt_database_path /default',
                    "bash",
                    str(CLI),
                ],
                check=True,
                capture_output=True,
                text=True,
                input="rel/db\n",
                cwd=tmp,
            )
            self.assertEqual(relative.stdout, f"{tmp}/rel/db\n")

    def test_config_writer_skips_write_when_content_is_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config_dir = Path(tmp)
            config = config_dir / "local.config"
            config.write_text("database_path=/db\ndatabase_profile=img\n")
            config_dir.chmod(0o555)
            try:
                unchanged = write_config(config, "/db", "img")
                changed = write_config(config, "/db", "curated")
            finally:
                config_dir.chmod(0o755)
        self.assertEqual(unchanged.returncode, 0)
        self.assertNotEqual(changed.returncode, 0)


class DatabaseInstallPolicyTests(unittest.TestCase):
    def ensure(
        self, validate_status: int, explicit_path: int = 0
    ) -> tuple[subprocess.CompletedProcess, list[str]]:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=validate_status)
            result = sandbox.run(f"ensure_database /database curated {explicit_path}")
            return result, sandbox.manager_calls()

    def installs(self, calls: list[str]) -> list[str]:
        return [call for call in calls if " install " in call]

    def test_validation_error_is_shown_to_the_user(self) -> None:
        for status in (2, 3, 4):
            with self.subTest(status=status):
                result, _ = self.ensure(status)
                self.assertIn(VALIDATE_MESSAGE, result.stderr)

    def test_absent_profile_is_installed_without_force(self) -> None:
        result, calls = self.ensure(3)
        installs = self.installs(calls)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(installs), 1)
        self.assertNotIn("--force", installs[0])

    def test_absent_profile_is_installed_under_explicit_database_path(self) -> None:
        result, calls = self.ensure(3, explicit_path=1)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.installs(calls)), 1)

    def test_damaged_profile_is_replaced_with_force(self) -> None:
        result, calls = self.ensure(4)
        installs = self.installs(calls)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(installs), 1)
        self.assertIn("--force", installs[0])

    def test_damaged_profile_under_explicit_database_path_is_not_replaced(self) -> None:
        result, calls = self.ensure(4, explicit_path=1)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.installs(calls), [])
        self.assertIn("--database_path", result.stderr)

    def test_other_validation_failure_stops_without_download(self) -> None:
        result, calls = self.ensure(2)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.installs(calls), [])


class SetupArgumentTests(unittest.TestCase):
    def test_help_prints_usage_without_side_effects(self) -> None:
        for flag in ("-h", "--help"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as tmp:
                sandbox = Sandbox(Path(tmp))
                result = sandbox.run("setup_database " + flag)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("Usage:", result.stdout)
                self.assertEqual(sandbox.manager_calls(), [])
                self.assertFalse(sandbox.config.exists())

    def test_invalid_arguments_fail_with_status_two_before_any_side_effect(
        self,
    ) -> None:
        cases = (
            ("--bogus",),
            ("extra",),
            ("--database_profile",),
            ("--database_profile", "--update"),
            ("--database_path",),
            ("--database_path=",),
            ("--database_path=-cache",),
            ("--database_profile=-curated",),
            ("--update=1",),
        )
        for arguments in cases:
            with (
                self.subTest(arguments=arguments),
                tempfile.TemporaryDirectory() as tmp,
            ):
                sandbox = Sandbox(Path(tmp))
                result = sandbox.run('setup_database "$@"', *arguments)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(sandbox.manager_calls(), [])
                self.assertFalse(sandbox.config.exists())

    def test_documented_argument_forms_are_accepted(self) -> None:
        cases = (
            ("--database_profile", "img"),
            ("--database_profile=img",),
            ("--database_path", "/db"),
            ("--database_path=/db",),
            ("--database_profile", "img", "--database_path", "/db", "--update"),
            (),
        )
        for arguments in cases:
            with (
                self.subTest(arguments=arguments),
                tempfile.TemporaryDirectory() as tmp,
            ):
                sandbox = Sandbox(Path(tmp), installed=True)
                result = sandbox.run('setup_database "$@"', *arguments)
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_unknown_profile_is_rejected_before_any_install(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=3)
            result = sandbox.run("setup_database --database_profile nope")
            calls = sandbox.manager_calls()
            self.assertFalse(sandbox.config.exists())
        self.assertEqual(result.returncode, 2)
        self.assertIn("nope", result.stderr)
        self.assertIn("curated", result.stderr)
        self.assertIn("img", result.stderr)
        self.assertEqual(len(calls), 1)
        self.assertIn(" profiles", calls[0])

    def test_config_is_written_after_the_database_is_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=3)
            result = sandbox.run("setup_database --database_profile img")
            config = sandbox.config.read_text()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(config, "database_path=/database\ndatabase_profile=img\n")

    def test_config_is_not_written_when_the_database_is_not_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=2)
            result = sandbox.run("setup_database --database_profile img")
            config_exists = sandbox.config.exists()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(config_exists)

    def test_update_fails_when_the_update_check_is_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(
                Path(tmp),
                installed=True,
                update_status="unavailable\t1.0.1\t-\toffline",
            )
            forced = sandbox.run("setup_database --database_profile img --update")
            plain = sandbox.run("setup_database --database_profile img")
        self.assertNotEqual(forced.returncode, 0)
        self.assertIn("offline", forced.stderr)
        self.assertEqual(plain.returncode, 0, plain.stderr)
        self.assertIn("offline", plain.stderr)

    def test_update_is_checked_after_a_fresh_install(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(
                Path(tmp),
                validate_status=3,
                update_status="unavailable\t1.0.1\t-\toffline",
            )
            forced = sandbox.run("setup_database --database_profile img --update")
            forced_calls = sandbox.manager_calls()
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=3)
            plain = sandbox.run("setup_database --database_profile img")
            plain_calls = sandbox.manager_calls()

        self.assertNotEqual(forced.returncode, 0)
        self.assertTrue(any("check-update" in call for call in forced_calls))
        self.assertEqual(plain.returncode, 0, plain.stderr)
        self.assertTrue(any("check-update" in call for call in plain_calls))

    def test_failed_config_write_stops_setup_and_keeps_the_saved_config(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), installed=True)
            sandbox.config.write_text("database_path=/old\ndatabase_profile=curated\n")
            result = sandbox.run(
                # The redirection still creates the temporary file.
                'printf() { [[ "$1" == database_path=* ]] && return 1; '
                'builtin printf "$@"; }; '
                "setup_database --database_profile img --database_path /new"
            )
            saved = sandbox.config.read_text()
            leftovers = [path.name for path in sandbox.config_dir.iterdir()]

        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(saved, "database_path=/old\ndatabase_profile=curated\n")
        self.assertEqual(leftovers, ["local.config"])


class RunArgumentTests(unittest.TestCase):
    def test_runs_use_quick_validation_and_cached_updates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), installed=True)
            result = sandbox.run("run_pipeline --query /q")
            calls = sandbox.manager_calls()
        self.assertEqual(result.returncode, 0, result.stderr)
        validation = next(call for call in calls if " validate " in call)
        update = next(call for call in calls if " check-update " in call)
        self.assertIn("--quick", validation)
        self.assertIn("--cached", update)

    def test_setup_always_uses_full_validation_and_fresh_update_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), installed=True)
            sandbox.env["SSUEXTRACT_SKIP_UPDATE_CHECK"] = "1"
            result = sandbox.run("setup_database --database_profile img")
            calls = sandbox.manager_calls()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn(
            "--quick", next(call for call in calls if " validate " in call)
        )
        self.assertNotIn(
            "--cached", next(call for call in calls if " check-update " in call)
        )

    def test_skip_update_check_keeps_run_and_example_ready_checks(self) -> None:
        for command in ("run_pipeline --query /q", "run_smoke"):
            with self.subTest(command=command), tempfile.TemporaryDirectory() as tmp:
                sandbox = Sandbox(Path(tmp), installed=True)
                sandbox.env["SSUEXTRACT_SKIP_UPDATE_CHECK"] = "1"
                result = sandbox.run(command)
                calls = sandbox.manager_calls()
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(any(" validate " in call for call in calls))
                self.assertFalse(any(" check-update " in call for call in calls))
                self.assertTrue(any(call.startswith("NEXTFLOW ") for call in calls))

    def test_example_stops_after_failed_preparation_or_version_lookup(self) -> None:
        for failed_step in ("mktemp", "cp", "version"):
            with (
                self.subTest(failed_step=failed_step),
                tempfile.TemporaryDirectory() as tmp,
            ):
                root = Path(tmp)
                sandbox = Sandbox(root, installed=True)
                stage = root / "stage"
                stage.mkdir()
                sandbox.env["SMOKE_STAGE"] = str(stage)
                prefix = ""
                if failed_step == "mktemp":
                    prefix = 'mktemp() { printf "%s\\n" "$SMOKE_STAGE"; return 71; }; '
                elif failed_step == "cp":
                    prefix = "cp() { return 71; }; "
                else:
                    sandbox.python.write_text(
                        FAKE_MANAGER.replace(
                            "version) echo 1.0.1",
                            'version) if grep -q "^NEXTFLOW " "$CALLS"; '
                            "then exit 71; fi; echo 1.0.1",
                        )
                    )
                result = sandbox.run(
                    prefix + "run_smoke --database_profile img || exit $?"
                )
                calls = sandbox.manager_calls()
                self.assertEqual(result.returncode, 71, result.stderr)
                self.assertFalse(
                    any("validate_example_output.py" in call for call in calls)
                )
                if failed_step != "version":
                    self.assertFalse(
                        any(call.startswith("NEXTFLOW ") for call in calls)
                    )

    def test_example_resolves_unsaved_path_once_for_run_and_version(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), installed=True)
            result = sandbox.run(
                'resolve_database_path() { printf "PATH_SELECTION\\n" >> "$CALLS"; '
                'printf "/selected database\\n"; }; run_smoke --database_profile img'
            )
            calls = sandbox.manager_calls()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(calls.count("PATH_SELECTION"), 1)
        self.assertTrue(
            any(
                "NEXTFLOW " in call and "--database_path /selected database" in call
                for call in calls
            )
        )
        versions = [call for call in calls if " version " in call]
        self.assertEqual(len(versions), 2)
        self.assertTrue(
            all("--root /selected database --profile img" in call for call in versions)
        )
        self.assertIn("--profile img --database-version 1.0.1", calls[-1])

    def test_database_options_without_a_value_fail_with_status_two(self) -> None:
        cases = (
            ("--database_profile",),
            ("--database_path",),
            ("--database_profile", "--query", "/queries"),
            ("--database_path", "--query", "/queries"),
            ("--database_path=",),
            ("--database_path=-cache",),
            ("--database_profile=-curated",),
        )
        for command in ("run_pipeline", "run_smoke"):
            for arguments in cases:
                with (
                    self.subTest(command=command, arguments=arguments),
                    tempfile.TemporaryDirectory() as tmp,
                ):
                    sandbox = Sandbox(Path(tmp))
                    result = sandbox.run(f'{command} "$@"', *arguments)
                    self.assertEqual(result.returncode, 2, result.stderr)
                    self.assertEqual(sandbox.manager_calls(), [])

    def test_unknown_profile_is_rejected_before_any_install(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=3)
            result = sandbox.run("run_pipeline --database_profile nope --query /q")
            calls = sandbox.manager_calls()
        self.assertEqual(result.returncode, 2)
        self.assertIn("curated", result.stderr)
        self.assertEqual(len(calls), 1)
        self.assertIn(" profiles", calls[0])

    def test_run_writes_config_only_after_the_database_is_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=3)
            ok = sandbox.run("run_pipeline --query /q")
            config = sandbox.config.read_text()
            nextflow_calls = [
                c for c in sandbox.manager_calls() if c.startswith("NEXTFLOW")
            ]
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertEqual(config, "database_path=/database\ndatabase_profile=curated\n")
        self.assertEqual(len(nextflow_calls), 1)

    def test_run_does_not_write_config_when_the_database_is_not_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=2)
            result = sandbox.run("run_pipeline --query /q")
            config_exists = sandbox.config.exists()
            nextflow_calls = [
                c for c in sandbox.manager_calls() if c.startswith("NEXTFLOW")
            ]
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(config_exists)
        self.assertEqual(nextflow_calls, [])

    def test_run_with_explicit_database_path_does_not_replace_a_damaged_profile(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sandbox = Sandbox(Path(tmp), validate_status=4)
            result = sandbox.run(
                "run_pipeline --database_path /custom --database_profile img --query /q"
            )
            calls = sandbox.manager_calls()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(VALIDATE_MESSAGE, result.stderr)
        self.assertFalse(any(" install " in call for call in calls))


if __name__ == "__main__":
    unittest.main()

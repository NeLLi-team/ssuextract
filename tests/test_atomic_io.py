import errno
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import atomic_io


class DurableReplaceTests(unittest.TestCase):
    def test_replace_syncs_parent_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "staged"
            destination = root / "release"
            source.write_text("complete", encoding="utf-8")
            with mock.patch.object(atomic_io.os, "fsync") as fsync:
                atomic_io.replace_and_fsync(source, destination)

            self.assertEqual(fsync.call_count, 2)
            self.assertFalse(source.exists())
            self.assertEqual(destination.read_text(encoding="utf-8"), "complete")

    def test_unsupported_directory_sync_keeps_atomic_replace_portable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "staged"
            destination = root / "release"
            source.write_text("complete", encoding="utf-8")
            unsupported = OSError(errno.EINVAL, "directory fsync unsupported")

            with mock.patch.object(atomic_io, "fsync_file"), mock.patch.object(
                atomic_io.os, "fsync", side_effect=unsupported
            ):
                atomic_io.replace_and_fsync(source, destination)

            self.assertEqual(destination.read_text(encoding="utf-8"), "complete")

    def test_directory_replace_syncs_files_and_directories(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "staged"
            nested = source / "tables"
            nested.mkdir(parents=True)
            (nested / "taxonomy.parquet").write_text("complete", encoding="utf-8")
            destination = root / "release"

            with mock.patch.object(atomic_io, "fsync_file") as fsync_file, mock.patch.object(
                atomic_io, "fsync_directory"
            ) as fsync_directory:
                atomic_io.replace_and_fsync(source, destination)

            fsync_file.assert_called_once()
            self.assertGreaterEqual(fsync_directory.call_count, 3)
            self.assertTrue((destination / "tables" / "taxonomy.parquet").is_file())


if __name__ == "__main__":
    unittest.main()

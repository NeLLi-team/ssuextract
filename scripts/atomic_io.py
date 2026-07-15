"""Durable installation primitive for atomically staged files."""

from __future__ import annotations

import errno
import os
from pathlib import Path


def fsync_file(path: str | Path) -> None:
    """Persist one regular file before publishing its directory entry."""

    file_fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(file_fd)
    finally:
        os.close(file_fd)


def fsync_directory(path: str | Path) -> None:
    """Persist directory-entry changes where the filesystem supports it."""

    directory_fd = os.open(path, os.O_RDONLY)
    try:
        try:
            os.fsync(directory_fd)
        except OSError as error:
            if error.errno not in {errno.EINVAL, errno.ENOTSUP}:
                raise
    finally:
        os.close(directory_fd)


def replace_and_fsync(source: str | Path, destination: str | Path) -> None:
    """Replace a staged file and sync its directory entry before returning."""

    staged = Path(source)
    target = Path(destination)
    source_parent = staged.parent
    if staged.is_dir():
        for file_path in sorted(path for path in staged.rglob("*") if path.is_file()):
            fsync_file(file_path)
        for directory in sorted(
            (path for path in staged.rglob("*") if path.is_dir()),
            key=lambda path: len(path.parts),
            reverse=True,
        ):
            fsync_directory(directory)
        fsync_directory(staged)
    else:
        fsync_file(staged)
    os.replace(staged, target)
    fsync_directory(target.parent)
    if source_parent != target.parent:
        fsync_directory(source_parent)

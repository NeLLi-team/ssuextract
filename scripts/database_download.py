"""Resumable, checksum-verified HTTP download for database archives."""

from __future__ import annotations

import hashlib
import http.client
import os
import re
import socket
import ssl
import stat
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import BinaryIO, Callable
from urllib.parse import urlparse

from atomic_io import fsync_directory


CHUNK_SIZE = 1024 * 1024
NETWORK_READ_SIZE = 64 * 1024
DEFAULT_ATTEMPTS = 5
DEFAULT_TIMEOUT = 60.0
DEFAULT_PROGRESS_INTERVAL = 1.0
_CONTENT_RANGE = re.compile(r"^bytes ([0-9]+)-([0-9]+)/([0-9]+)$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_RETRYABLE_HTTP_STATUS = {408, 429, 500, 502, 503, 504}

ProgressReporter = Callable[[str], None] | None


class DownloadError(RuntimeError):
    """The archive could not be transferred safely."""


class DownloadIntegrityError(DownloadError):
    """The response or completed archive violated its integrity contract."""


def _human_size(size: int) -> str:
    value = float(size)
    units = ("B", "KiB", "MiB", "GiB")
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    raise AssertionError("unreachable")


def _report(progress: ProgressReporter, message: str) -> None:
    if progress is not None:
        progress(message)


def _format_eta(seconds: float) -> str:
    rounded = max(0, int(seconds + 0.5))
    hours, remainder = divmod(rounded, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{seconds:02d}"
    return f"{minutes:02d}:{seconds:02d}"


def _progress_message(size: int, expected_size: int, rate: float) -> str:
    width = 24
    filled = min(width, size * width // expected_size)
    bar = "=" * filled + "-" * (width - filled)
    percent = size * 100 / expected_size
    if rate > 0:
        rate_text = f"{_human_size(int(rate))}/s"
        eta_text = _format_eta((expected_size - size) / rate)
    else:
        rate_text = "-- B/s"
        eta_text = "--:--"
    return (
        f"Downloading archive: [{bar}] {percent:5.1f}% "
        f"({_human_size(size)} of {_human_size(expected_size)}) "
        f"{rate_text} ETA {eta_text}"
    )


def _read_available(response: BinaryIO) -> bytes:
    read1 = getattr(response, "read1", None)
    if callable(read1):
        return read1(NETWORK_READ_SIZE)
    return response.read(NETWORK_READ_SIZE)


def _header(response: BinaryIO, name: str) -> str | None:
    headers = getattr(response, "headers", None)
    if headers is None:
        return None
    return headers.get(name)


def _status(response: BinaryIO) -> int:
    status = getattr(response, "status", None)
    return 200 if status is None else int(status)


def _validate_response(
    response: BinaryIO,
    request_url: str,
    offset: int,
    expected_size: int,
) -> tuple[bool, int]:
    """Validate response framing and return restart state plus response bytes."""

    final_url = response.geturl() if hasattr(response, "geturl") else request_url
    parsed = urlparse(final_url)
    if parsed.scheme != "https" or not parsed.netloc:
        raise DownloadIntegrityError(
            f"Database archive redirected to a non-HTTPS URL: {final_url}"
        )

    status = _status(response)
    content_range = _header(response, "Content-Range")
    content_length = _header(response, "Content-Length")
    content_encoding = _header(response, "Content-Encoding")
    if content_encoding is not None and content_encoding.lower() != "identity":
        raise DownloadIntegrityError(
            f"Archive response uses unsupported Content-Encoding: {content_encoding}"
        )
    if offset == 0:
        if status != 200:
            raise DownloadIntegrityError(
                f"Initial archive request returned HTTP {status}; expected HTTP 200"
            )
        if content_range is not None:
            raise DownloadIntegrityError(
                "Initial archive response unexpectedly included Content-Range"
            )
        expected_response_size = expected_size
        ignored_range = False
    elif status == 200:
        if content_range is not None:
            raise DownloadIntegrityError(
                "Archive server returned HTTP 200 with an unexpected Content-Range"
            )
        expected_response_size = expected_size
        ignored_range = True
    elif status == 206:
        match = _CONTENT_RANGE.fullmatch(content_range or "")
        if match is None:
            raise DownloadIntegrityError(
                "Archive range response is missing a valid Content-Range"
            )
        start, end, total = (int(value) for value in match.groups())
        if start != offset or end < start or end >= total or total != expected_size:
            raise DownloadIntegrityError(
                f"Archive range response does not match bytes {offset}-/{expected_size}"
            )
        expected_response_size = end - start + 1
        ignored_range = False
    else:
        raise DownloadIntegrityError(
            f"Archive range request returned HTTP {status}; expected HTTP 206"
        )

    if content_length is not None:
        try:
            actual_response_size = int(content_length)
        except ValueError as error:
            raise DownloadIntegrityError(
                f"Archive response has invalid Content-Length: {content_length!r}"
            ) from error
        if actual_response_size != expected_response_size:
            raise DownloadIntegrityError(
                f"Archive response length mismatch: expected {expected_response_size} bytes, "
                f"server reported {actual_response_size}"
            )
    return ignored_range, expected_response_size


def _retryable(error: BaseException) -> bool:
    if isinstance(error, ssl.SSLCertVerificationError):
        return False
    if isinstance(error, urllib.error.HTTPError):
        return error.code in _RETRYABLE_HTTP_STATUS
    if isinstance(error, urllib.error.URLError):
        reason = error.reason
        if isinstance(reason, ssl.SSLCertVerificationError):
            return False
        return True
    return isinstance(
        error,
        (
            TimeoutError,
            socket.timeout,
            ConnectionError,
            http.client.IncompleteRead,
            http.client.RemoteDisconnected,
            OSError,
        ),
    )


def _one_line(error: BaseException) -> str:
    return " ".join(str(error).split()) or error.__class__.__name__


def _hash_partial(path: Path) -> tuple[int, "hashlib._Hash"]:
    digest = hashlib.sha256()
    size = 0
    if path.exists():
        if not path.is_file() or path.is_symlink():
            raise DownloadError(f"Partial archive is not a regular file: {path}")
        try:
            with path.open("rb") as handle:
                while chunk := handle.read(CHUNK_SIZE):
                    size += len(chunk)
                    digest.update(chunk)
        except OSError as error:
            raise DownloadError(f"Could not read partial archive {path}: {error}") from error
    return size, digest


def _remove_partial(path: Path) -> None:
    try:
        path.unlink(missing_ok=True)
        fsync_directory(path.parent)
    except OSError as error:
        raise DownloadError(f"Could not remove invalid partial archive {path}: {error}") from error


def _open_output(path: Path, size: int):
    flags = os.O_WRONLY | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0)
    flags |= os.O_APPEND if size else os.O_TRUNC
    try:
        descriptor = os.open(path, flags, 0o600)
    except OSError as error:
        raise DownloadError(f"Could not open partial archive {path}: {error}") from error
    try:
        details = os.fstat(descriptor)
        if not stat.S_ISREG(details.st_mode):
            raise DownloadError(f"Partial archive is not a regular file: {path}")
        if details.st_size != size:
            raise DownloadError(
                f"Partial archive changed during download: expected {size} bytes, "
                f"found {details.st_size}"
            )
        return os.fdopen(descriptor, "ab" if size else "wb")
    except BaseException:
        os.close(descriptor)
        raise


def download_verified_archive(
    url: str,
    destination: str | Path,
    expected_size: int,
    expected_sha256: str,
    opener: Callable[..., BinaryIO] = urllib.request.urlopen,
    progress: ProgressReporter = None,
    timeout: float = DEFAULT_TIMEOUT,
    attempts: int = DEFAULT_ATTEMPTS,
    sleeper: Callable[[float], None] | None = None,
    clock: Callable[[], float] | None = None,
    progress_interval: float = DEFAULT_PROGRESS_INTERVAL,
) -> None:
    """Download an HTTPS archive, retaining and resuming a release-bound partial file."""

    parsed = urlparse(url)
    if parsed.scheme != "https" or not parsed.netloc:
        raise ValueError("database archive URL must use HTTPS")
    if expected_size <= 0 or _SHA256.fullmatch(expected_sha256) is None:
        raise ValueError("expected archive size and SHA-256 are invalid")
    if timeout <= 0 or attempts <= 0 or progress_interval <= 0:
        raise ValueError("download timeout, attempts, and progress interval must be positive")
    destination_path = Path(destination)
    try:
        destination_path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise DownloadError(
            f"Could not create download directory {destination_path.parent}: {error}"
        ) from error
    size, digest = _hash_partial(destination_path)
    if size > expected_size:
        _report(progress, "Discarding an oversized partial database archive.")
        _remove_partial(destination_path)
        size, digest = 0, hashlib.sha256()
    elif size == expected_size:
        if digest.hexdigest() == expected_sha256:
            _report(progress, f"Using verified downloaded archive at {destination_path}.")
            return
        _report(progress, "Discarding a complete partial archive with the wrong SHA-256.")
        _remove_partial(destination_path)
        size, digest = 0, hashlib.sha256()

    resumed_existing = size > 0
    clean_restart_used = False
    range_restart_used = False
    sleep = time.sleep if sleeper is None else sleeper
    monotonic = time.monotonic if clock is None else clock
    attempt = 1
    while True:
        percent = size * 100 // expected_size
        if size:
            _report(
                progress,
                f"Resuming archive at {_human_size(size)} of {_human_size(expected_size)} "
                f"({percent}%).",
            )
        else:
            _report(progress, _progress_message(0, expected_size, 0))

        request = urllib.request.Request(url)
        if size:
            request.add_header("Range", f"bytes={size}-{expected_size - 1}")

        try:
            response = opener(request, timeout=timeout)
        except (OSError, ValueError, http.client.HTTPException) as error:
            status = error.code if isinstance(error, urllib.error.HTTPError) else None
            if isinstance(error, urllib.error.HTTPError):
                error.close()
            if status == 416 and size and not range_restart_used:
                _report(
                    progress,
                    "Archive server rejected the retained byte range; discarding the "
                    "partial archive and restarting once from byte 0.",
                )
                _remove_partial(destination_path)
                size, digest = 0, hashlib.sha256()
                range_restart_used = True
                resumed_existing = False
                attempt = 1
                continue
            if _retryable(error) and attempt < attempts:
                delay = min(2 ** (attempt - 1), 16)
                _report(
                    progress,
                    f"Download interrupted ({_one_line(error)}); retrying in {delay}s "
                    f"from {_human_size(size)} (attempt {attempt + 1} of {attempts}).",
                )
                sleep(delay)
                attempt += 1
                continue
            detail = _one_line(error)
            raise DownloadError(
                f"Could not download database archive {url}: {detail}. Partial download "
                f"retained at {destination_path} ({size} of {expected_size} bytes); "
                "rerun setup to resume."
            ) from error

        network_error: BaseException | None = None
        try:
            with response:
                ignored_range, expected_response_bytes = _validate_response(
                    response, url, size, expected_size
                )
                if ignored_range:
                    _report(
                        progress,
                        "Archive server ignored the byte range; restarting safely from byte 0.",
                    )
                    size, digest = 0, hashlib.sha256()
                progress_started = monotonic()
                progress_start_size = size
                last_progress = progress_started
                if size or ignored_range:
                    _report(progress, _progress_message(size, expected_size, 0))
                response_bytes = 0
                with _open_output(destination_path, size) as output:
                    while True:
                        try:
                            chunk = _read_available(response)
                        except (
                            OSError,
                            http.client.HTTPException,
                        ) as error:
                            if _retryable(error):
                                network_error = error
                                break
                            raise DownloadError(
                                f"Could not read database archive {url}: {_one_line(error)}"
                            ) from error
                        if not chunk:
                            break
                        if response_bytes + len(chunk) > expected_response_bytes:
                            raise DownloadIntegrityError(
                                "Archive response exceeded its declared byte range"
                            )
                        if size + len(chunk) > expected_size:
                            raise DownloadIntegrityError(
                                f"Archive size exceeds catalog value: expected "
                                f"{expected_size} bytes"
                            )
                        try:
                            output.write(chunk)
                        except OSError as error:
                            raise DownloadError(
                                f"Could not write partial archive {destination_path}: {error}"
                            ) from error
                        size += len(chunk)
                        response_bytes += len(chunk)
                        digest.update(chunk)
                        now = monotonic()
                        elapsed = now - progress_started
                        if size == expected_size or now - last_progress >= progress_interval:
                            rate = (
                                (size - progress_start_size) / elapsed
                                if elapsed > 0
                                else 0
                            )
                            _report(
                                progress,
                                _progress_message(size, expected_size, rate),
                            )
                            last_progress = now
                    output.flush()
                    os.fsync(output.fileno())
                fsync_directory(destination_path.parent)
        except DownloadError:
            raise
        except OSError as error:
            raise DownloadError(
                f"Could not persist partial archive {destination_path}: {error}"
            ) from error

        if network_error is None and response_bytes != expected_response_bytes:
            network_error = http.client.IncompleteRead(
                b"", expected_response_bytes - response_bytes
            )

        if network_error is None and size < expected_size:
            _report(
                progress,
                f"Archive server completed a sub-range; continuing from "
                f"{_human_size(size)}.",
            )
            attempt = 1
            continue

        if network_error is None and size == expected_size:
            actual_sha256 = digest.hexdigest()
            if actual_sha256 == expected_sha256:
                return
            if resumed_existing and not clean_restart_used:
                _report(
                    progress,
                    "Resumed archive failed SHA-256 verification; discarding it and "
                    "retrying once from byte 0.",
                )
                _remove_partial(destination_path)
                size, digest = 0, hashlib.sha256()
                clean_restart_used = True
                resumed_existing = False
                attempt = 1
                continue
            _remove_partial(destination_path)
            raise DownloadIntegrityError(
                f"Archive SHA-256 mismatch: expected {expected_sha256}, "
                f"found {actual_sha256}"
            )

        if attempt >= attempts:
            reason = (
                _one_line(network_error)
                if network_error is not None
                else "the server ended the response before the expected size"
            )
            raise DownloadError(
                f"Database archive download stopped after {attempts} attempts ({reason}). "
                f"Partial download retained at {destination_path} "
                f"({size} of {expected_size} bytes); rerun setup to resume."
            )
        delay = min(2 ** (attempt - 1), 16)
        reason = (
            _one_line(network_error)
            if network_error is not None
            else "server ended the response early"
        )
        _report(
            progress,
            f"Download interrupted ({reason}); retrying in {delay}s from "
            f"{_human_size(size)} (attempt {attempt + 1} of {attempts}).",
        )
        sleep(delay)
        attempt += 1

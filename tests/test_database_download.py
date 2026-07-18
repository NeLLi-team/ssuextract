import io
import ssl
import sys
import tempfile
import unittest
import urllib.error
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import database_download as downloader


class FakeResponse(io.BytesIO):
    def __init__(
        self,
        data: bytes,
        *,
        status: int = 200,
        headers: dict[str, str] | None = None,
        url: str = "https://example.org/database.tar.zst",
    ) -> None:
        super().__init__(data)
        self.status = status
        self.headers = headers or {}
        self.url = url

    def geturl(self) -> str:
        return self.url

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class InterruptingResponse(FakeResponse):
    def __init__(self, data: bytes, interrupt_after: int) -> None:
        super().__init__(data, headers={"Content-Length": str(len(data))})
        self.interrupt_after = interrupt_after
        self.interrupted = False

    def read(self, size: int = -1) -> bytes:
        if self.tell() >= self.interrupt_after and not self.interrupted:
            self.interrupted = True
            raise TimeoutError("simulated read timeout")
        remaining = self.interrupt_after - self.tell()
        return super().read(min(size, remaining))

    def read1(self, size: int = -1) -> bytes:
        return self.read(size)


class AvailableBytesResponse(FakeResponse):
    def __init__(self, data: bytes, available: int) -> None:
        super().__init__(data, headers={"Content-Length": str(len(data))})
        self.available = available
        self.read1_calls = 0

    def read(self, size: int = -1) -> bytes:
        raise TimeoutError("full-buffer read waited too long")

    def read1(self, size: int = -1) -> bytes:
        self.read1_calls += 1
        return super().read(min(size, self.available))


def sha256(data: bytes) -> str:
    import hashlib

    return hashlib.sha256(data).hexdigest()


def range_response(data: bytes, start: int, end: int | None = None) -> FakeResponse:
    end = len(data) - 1 if end is None else end
    return FakeResponse(
        data[start : end + 1],
        status=206,
        headers={
            "Content-Range": f"bytes {start}-{end}/{len(data)}",
            "Content-Length": str(end - start + 1),
        },
    )


class DatabaseDownloadTests(unittest.TestCase):
    def test_input_and_redirect_urls_must_remain_https(self) -> None:
        data = b"abcdefghij"
        with tempfile.TemporaryDirectory() as tmp, self.assertRaisesRegex(
            ValueError, "must use HTTPS"
        ):
            downloader.download_verified_archive(
                "http://example.org/database.tar.zst",
                Path(tmp) / "database.part",
                len(data),
                sha256(data),
            )

        with tempfile.TemporaryDirectory() as tmp, self.assertRaisesRegex(
            downloader.DownloadIntegrityError, "non-HTTPS"
        ):
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                Path(tmp) / "database.part",
                len(data),
                sha256(data),
                opener=lambda request, timeout: FakeResponse(
                    data,
                    headers={"Content-Length": str(len(data))},
                    url="http://mirror.example.org/database.tar.zst",
                ),
            )

    def test_certificate_failure_is_not_retried(self) -> None:
        data = b"abcdefghij"
        calls = 0

        def opener(request, timeout):
            nonlocal calls
            calls += 1
            raise urllib.error.URLError(
                ssl.SSLCertVerificationError("simulated certificate failure")
            )

        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(
            downloader.DownloadError
        ):
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                Path(tmp) / "database.part",
                len(data),
                sha256(data),
                opener=opener,
                sleeper=lambda delay: self.fail("certificate failure must not retry"),
            )
        self.assertEqual(calls, 1)

    def test_encoded_response_is_rejected(self) -> None:
        data = b"abcdefghij"
        with tempfile.TemporaryDirectory() as tmp, self.assertRaisesRegex(
            downloader.DownloadIntegrityError, "Content-Encoding"
        ):
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                Path(tmp) / "database.part",
                len(data),
                sha256(data),
                opener=lambda request, timeout: FakeResponse(
                    data,
                    headers={
                        "Content-Length": str(len(data)),
                        "Content-Encoding": "gzip",
                    },
                ),
            )

    def test_midstream_timeout_resumes_within_one_invocation(self) -> None:
        data = b"0123456789"
        requests = []

        def opener(request, timeout):
            requests.append(request)
            if len(requests) == 1:
                return InterruptingResponse(data, 4)
            return range_response(data, 4)

        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                destination,
                len(data),
                sha256(data),
                opener=opener,
                sleeper=lambda delay: None,
            )
            self.assertEqual(destination.read_bytes(), data)

        self.assertIsNone(requests[0].get_header("Range"))
        self.assertEqual(requests[1].get_header("Range"), "bytes=4-9")

    def test_available_bytes_are_persisted_without_full_buffer_read(self) -> None:
        data = b"0123456789"
        response = AvailableBytesResponse(data, available=2)
        messages: list[str] = []
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                destination,
                len(data),
                sha256(data),
                opener=lambda request, timeout: response,
                progress=messages.append,
            )
            self.assertEqual(destination.read_bytes(), data)
        self.assertGreater(response.read1_calls, 1)
        self.assertTrue(any("100.0%" in message for message in messages))
        self.assertTrue(any("ETA" in message for message in messages))

    def test_progress_clock_controls_rate_eta_and_update_interval(self) -> None:
        data = b"012345"
        response = AvailableBytesResponse(data, available=2)
        messages: list[str] = []
        ticks = iter((0.0, 1.0, 2.0, 3.0))
        with tempfile.TemporaryDirectory() as tmp:
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                Path(tmp) / "database.part",
                len(data),
                sha256(data),
                opener=lambda request, timeout: response,
                progress=messages.append,
                clock=lambda: next(ticks),
                progress_interval=1.0,
            )
        self.assertTrue(
            any(
                "33.3%" in message and "2 B/s ETA 00:02" in message
                for message in messages
            )
        )
        self.assertEqual(sum("2 B/s" in message for message in messages), 3)

    def test_later_invocation_resumes_retained_partial(self) -> None:
        data = b"abcdefghij"
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            with self.assertRaisesRegex(downloader.DownloadError, "rerun setup to resume"):
                downloader.download_verified_archive(
                    "https://example.org/database.tar.zst",
                    destination,
                    len(data),
                    sha256(data),
                    opener=lambda request, timeout: InterruptingResponse(data, 3),
                    attempts=1,
                )
            self.assertEqual(destination.read_bytes(), data[:3])

            requests = []

            def resume(request, timeout):
                requests.append(request)
                return range_response(data, 3)

            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                destination,
                len(data),
                sha256(data),
                opener=resume,
            )
            self.assertEqual(destination.read_bytes(), data)
            self.assertEqual(requests[0].get_header("Range"), "bytes=3-9")

    def test_ignored_range_restarts_without_appending_duplicate_bytes(self) -> None:
        data = b"abcdefghij"
        messages: list[str] = []
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            destination.write_bytes(data[:4])
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                destination,
                len(data),
                sha256(data),
                opener=lambda request, timeout: FakeResponse(
                    data, headers={"Content-Length": str(len(data))}
                ),
                progress=messages.append,
            )
            self.assertEqual(destination.read_bytes(), data)
        self.assertTrue(any("ignored the byte range" in message for message in messages))

    def test_malformed_range_response_fails_without_changing_partial(self) -> None:
        data = b"abcdefghij"
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            destination.write_bytes(data[:4])
            with self.assertRaisesRegex(
                downloader.DownloadIntegrityError, "valid Content-Range"
            ):
                downloader.download_verified_archive(
                    "https://example.org/database.tar.zst",
                    destination,
                    len(data),
                    sha256(data),
                    opener=lambda request, timeout: FakeResponse(
                        data[4:],
                        status=206,
                        headers={"Content-Length": "6"},
                    ),
                )
            self.assertEqual(destination.read_bytes(), data[:4])

    def test_range_response_cannot_exceed_declared_content_range(self) -> None:
        data = b"abcdefghij"
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            destination.write_bytes(data[:4])
            with self.assertRaisesRegex(
                downloader.DownloadIntegrityError, "declared byte range"
            ):
                downloader.download_verified_archive(
                    "https://example.org/database.tar.zst",
                    destination,
                    len(data),
                    sha256(data),
                    opener=lambda request, timeout: FakeResponse(
                        data[4:],
                        status=206,
                        headers={
                            "Content-Range": "bytes 4-5/10",
                            "Content-Length": "2",
                        },
                    ),
                )
            self.assertEqual(destination.read_bytes(), data[:4])

    def test_server_chosen_sub_ranges_do_not_consume_failure_budget(self) -> None:
        data = b"abcdefghij"
        requests = []
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            destination.write_bytes(data[:1])

            def opener(request, timeout):
                requests.append(request)
                start = int(request.get_header("Range").split("=")[1].split("-")[0])
                return range_response(data, start, start)

            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                destination,
                len(data),
                sha256(data),
                opener=opener,
                attempts=2,
                sleeper=lambda delay: self.fail("clean sub-ranges must not back off"),
            )
            self.assertEqual(destination.read_bytes(), data)
        self.assertEqual(len(requests), 9)

    def test_rejected_retained_range_gets_one_clean_restart(self) -> None:
        data = b"abcdefghij"
        requests = []
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            destination.write_bytes(data[:4])

            def opener(request, timeout):
                requests.append(request)
                if len(requests) == 1:
                    raise urllib.error.HTTPError(
                        request.full_url, 416, "range rejected", {}, None
                    )
                return FakeResponse(data, headers={"Content-Length": str(len(data))})

            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                destination,
                len(data),
                sha256(data),
                opener=opener,
            )
            self.assertEqual(destination.read_bytes(), data)
        self.assertEqual(requests[0].get_header("Range"), "bytes=4-9")
        self.assertIsNone(requests[1].get_header("Range"))

    def test_corrupt_retained_prefix_gets_one_clean_restart(self) -> None:
        data = b"abcdefghij"
        requests = []
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            destination.write_bytes(b"xxxx")

            def opener(request, timeout):
                requests.append(request)
                if len(requests) == 1:
                    return range_response(data, 4)
                return FakeResponse(data, headers={"Content-Length": str(len(data))})

            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                destination,
                len(data),
                sha256(data),
                opener=opener,
            )
            self.assertEqual(destination.read_bytes(), data)
        self.assertEqual(requests[0].get_header("Range"), "bytes=4-9")
        self.assertIsNone(requests[1].get_header("Range"))

    def test_complete_verified_partial_uses_no_network(self) -> None:
        data = b"abcdefghij"
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "database.part"
            destination.write_bytes(data)
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                destination,
                len(data),
                sha256(data),
                opener=lambda request, timeout: self.fail("network should not be used"),
            )

    def test_partial_archive_symlink_is_rejected(self) -> None:
        data = b"abcdefghij"
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "target"
            target.write_bytes(b"keep")
            destination = Path(tmp) / "database.part"
            destination.symlink_to(target)
            with self.assertRaisesRegex(
                downloader.DownloadError, "not a regular file"
            ):
                downloader.download_verified_archive(
                    "https://example.org/database.tar.zst",
                    destination,
                    len(data),
                    sha256(data),
                    opener=lambda request, timeout: FakeResponse(data),
                )
            self.assertEqual(target.read_bytes(), b"keep")

    def test_retryable_http_status_is_retried_but_404_is_not(self) -> None:
        data = b"abcdefghij"
        calls = 0

        def transient(request, timeout):
            nonlocal calls
            calls += 1
            if calls == 1:
                raise urllib.error.HTTPError(
                    request.full_url, 503, "unavailable", {}, None
                )
            return FakeResponse(data, headers={"Content-Length": str(len(data))})

        with tempfile.TemporaryDirectory() as tmp:
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                Path(tmp) / "database.part",
                len(data),
                sha256(data),
                opener=transient,
                sleeper=lambda delay: None,
            )
        self.assertEqual(calls, 2)

        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(
            downloader.DownloadError
        ):
            downloader.download_verified_archive(
                "https://example.org/database.tar.zst",
                Path(tmp) / "database.part",
                len(data),
                sha256(data),
                opener=lambda request, timeout: (_ for _ in ()).throw(
                    urllib.error.HTTPError(request.full_url, 404, "missing", {}, None)
                ),
                sleeper=lambda delay: self.fail("404 must not be retried"),
            )


if __name__ == "__main__":
    unittest.main()

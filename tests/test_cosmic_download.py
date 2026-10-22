"""Check archive integrity, safe publication and authentication separation offline."""
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import urllib.error
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "01_data_downloading_and_preprocessing"))
import download_cosmic_variants as cosmic


class CosmicDownloadTests(unittest.TestCase):
    def setUp(self):
        self.payload = b"test archive payload"
        self.record = {"filename": "Cosmic_Test.tar", "s3_bucket": "cosmic-downloads-production",
                       "api_bucket": "downloads",
                       "s3_object": "grch38/cosmic/v101/Cosmic_Test.tar",
                       "sha256": hashlib.sha256(self.payload).hexdigest()}
        self.url = "https://cosmic-downloads-production.s3.eu-west-2.amazonaws.com/" + self.record["s3_object"] + "?signature=SECRET"

    def test_api_auth_is_separate_from_storage_request(self):
        opener = Mock()
        opener.open.return_value = io.BytesIO(json.dumps({"url": self.url}).encode())
        url = cosmic.obtain_url(opener, "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted", self.record, "Basic PRIVATE")
        request = opener.open.call_args.args[0]
        self.assertEqual(request.get_header("Authorization"), "Basic PRIVATE")
        self.assertIn("bucket=downloads", request.full_url)
        opener.open.return_value = io.BytesIO(self.payload)
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / self.record["filename"]
            cosmic.download(opener, url, self.record, target)
            self.assertEqual(target.read_bytes(), self.payload)
            self.assertIsNone(opener.open.call_args.args[0].get_header("Authorization"))

    def test_bad_payload_never_becomes_an_archive(self):
        opener = Mock()
        opener.open.return_value = io.BytesIO(b"<html>login required</html>")
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / self.record["filename"]
            with self.assertRaisesRegex(cosmic.DownloadError, "SHA256 mismatch"):
                cosmic.download(opener, self.url, self.record, target)
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_existing_output_is_never_replaced(self):
        opener = Mock()
        opener.open.return_value = io.BytesIO(self.payload)
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / self.record["filename"]
            target.write_bytes(b"existing research data")
            with self.assertRaises(cosmic.DownloadError):
                cosmic.download(opener, self.url, self.record, target)
            self.assertEqual(target.read_bytes(), b"existing research data")
            self.assertEqual(len(list(Path(directory).iterdir())), 1)

    def test_http_errors_do_not_disclose_signed_urls(self):
        opener = Mock()
        opener.open.side_effect = urllib.error.HTTPError(self.url, 403, "Forbidden", {}, None)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(cosmic.DownloadError) as caught:
                cosmic.download(opener, self.url, self.record, Path(directory) / self.record["filename"])
            self.assertNotIn("SECRET", str(caught.exception))
            self.assertNotIn(self.url, str(caught.exception))

    def test_api_rejects_login_html_and_unknown_response(self):
        for body in [b"<html>login</html>", b"[]", b'{"message":"login required"}']:
            opener = Mock()
            opener.open.return_value = io.BytesIO(body)
            with self.subTest(body=body), self.assertRaises(cosmic.DownloadError):
                cosmic.obtain_url(opener, "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted", self.record, "Basic PRIVATE")

    def test_rejects_unrelated_hosts_plaintext_and_wrong_filename(self):
        for url in ["http://cog.sanger.ac.uk/Cosmic_Test.tar", "https://example.com/Cosmic_Test.tar",
                    "https://cosmic-downloads-production.s3.eu-west-2.amazonaws.com.evil.org/Cosmic_Test.tar",
                    "https://user:password@cog.sanger.ac.uk/Cosmic_Test.tar",
                    "https://cog.sanger.ac.uk/wrong.tar"]:
            with self.subTest(url=url), self.assertRaises(cosmic.DownloadError):
                cosmic.checked_url(url, self.record)

    def test_redirect_never_forwards_authorization(self):
        request = cosmic.urllib.request.Request("https://cancer.sanger.ac.uk/api", headers={"Authorization": "Basic PRIVATE"})
        with self.assertRaises(cosmic.DownloadError):
            cosmic.NoRedirect().redirect_request(request, None, 302, "Found", {}, "https://other.example")

    def test_verify_only_is_offline_and_receipt_contains_no_auth(self):
        manifest = {"release": "v101", "genome": "GRCh38", "release_selection": {"status": "candidate"}, "files": [self.record]}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            source.write_text(json.dumps(manifest))
            out = root / "raw"
            out.mkdir()
            (out / self.record["filename"]).write_bytes(self.payload)
            with patch.object(cosmic, "MANIFEST", source), patch.object(cosmic.urllib.request, "build_opener") as network:
                self.assertEqual(cosmic.main(["--out-dir", str(out), "--verify-only"]), 0)
                network.return_value.open.assert_not_called()
            receipt = json.loads(next(out.glob("download_verification_*.json")).read_text())
            self.assertTrue(receipt["download_verification_passed"])
            self.assertFalse(receipt["historical_release_confirmed"])


if __name__ == "__main__":
    unittest.main()

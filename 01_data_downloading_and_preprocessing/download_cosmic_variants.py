#!/usr/bin/env python3
"""Download pinned COSMIC archives with interactive authentication and SHA256 checks."""
import argparse
import base64
from datetime import datetime, timezone
import getpass
import hashlib
import json
import os
from pathlib import Path
import sys
import urllib.error
import urllib.parse
import urllib.request
import uuid


MANIFEST = Path(__file__).with_name("cosmic_v101_grch38.json")


class DownloadError(Exception):
    """An actionable error that contains no credentials or signed URL."""


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        # In particular, never forward the Basic Authorization header elsewhere.
        raise DownloadError("Unexpected redirect; check the official download instructions.")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_url(url, record):
    if not isinstance(url, str):
        raise DownloadError("The API did not return a download URL.")
    parsed = urllib.parse.urlsplit(url)
    host = parsed.hostname or ""
    bucket = record["s3_bucket"]
    trusted = host == "cog.sanger.ac.uk" or host == bucket + ".s3.amazonaws.com"
    trusted = trusted or (host.startswith(bucket + ".s3.") and host.endswith(".amazonaws.com"))
    if parsed.scheme != "https" or parsed.username or parsed.password or not trusted:
        raise DownloadError("Download URL is not an expected COSMIC HTTPS storage address.")
    if parsed.port not in (None, 443):
        raise DownloadError("Unexpected download URL port.")
    if Path(urllib.parse.unquote(parsed.path)).name != record["filename"]:
        raise DownloadError("Download URL does not identify the selected archive.")
    return url


def obtain_url(opener, endpoint, record, authorization):
    # The public API bucket alias differs from the physical storage bucket name.
    query = urllib.parse.urlencode({"path": record["s3_object"], "bucket": record["api_bucket"]})
    request = urllib.request.Request(endpoint + "?" + query, headers={
        "Authorization": authorization, "Accept": "application/json"})
    try:
        with opener.open(request, timeout=60) as response:
            body = response.read(1024 * 1024)
        payload = json.loads(body)
    except urllib.error.HTTPError as error:
        raise DownloadError("COSMIC API returned HTTP %s. Check account login/download access "
                            "and the current Scripted download instructions." % error.code) from None
    except (urllib.error.URLError, TimeoutError, OSError, ValueError):
        raise DownloadError("Could not read a JSON download link from the COSMIC API.") from None
    if not isinstance(payload, dict):
        raise DownloadError("Unexpected COSMIC API response; expected a JSON object with 'url'.")
    return checked_url(payload.get("url"), record)


def download(opener, url, record, destination):
    checked_url(url, record)
    temporary = destination.with_name(destination.name + "." + uuid.uuid4().hex + ".part")
    try:
        # This request deliberately has no Authorization header.
        request = urllib.request.Request(url)
        with opener.open(request, timeout=60) as source, temporary.open("xb") as target:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                target.write(block)
        if sha256(temporary) != record["sha256"]:
            raise DownloadError("SHA256 mismatch for %s; no final archive was saved." % record["filename"])
        # Hard-link publication refuses an existing destination, including concurrent writers.
        os.link(temporary, destination)
    except urllib.error.HTTPError as error:
        raise DownloadError("Storage returned HTTP %s; obtain a fresh download link." % error.code) from None
    except (urllib.error.URLError, TimeoutError, OSError):
        raise DownloadError("Download or file publication failed for %s; check network, "
                            "disk space and output path." % record["filename"]) from None
    finally:
        temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, help="Raw input directory outside the Git checkout")
    parser.add_argument("--plan", action="store_true", help="Print pinned public metadata; no login or network")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--verify-only", action="store_true", help="Check archives already downloaded in a browser")
    modes.add_argument("--signed-links", action="store_true", help="Paste temporary URLs at hidden prompts instead of account login")
    args = parser.parse_args(argv)
    manifest = json.loads(MANIFEST.read_text())
    if args.plan:
        print(json.dumps(manifest, indent=2))
        return 0
    if args.out_dir is None:
        parser.error("--out-dir is required unless --plan is used")
    pending = []
    for record in manifest["files"]:
        destination = args.out_dir / record["filename"]
        if destination.exists():
            if sha256(destination) != record["sha256"]:
                raise DownloadError("Existing archive has the wrong SHA256: %s. Preserve it "
                                    "and choose a new output directory." % record["filename"])
            print(record["filename"] + ": SHA256 OK (existing)")
        else:
            pending.append(record)
    if args.verify_only and pending:
        raise DownloadError("Missing archives: " + ", ".join(r["filename"] for r in pending))
    authorization = None
    if pending and not args.signed_links:
        if not sys.stdin.isatty():
            raise DownloadError("Run in an interactive terminal for hidden account input, "
                                "or download in a browser and use --verify-only.")
        email = input("COSMIC account email: ").strip()
        password = getpass.getpass("COSMIC password (not saved): ")
        if not email or not password:
            raise DownloadError("Both account email and password are required.")
        authorization = "Basic " + base64.b64encode((email + ":" + password).encode()).decode("ascii")
        del password
    if pending and args.signed_links and not sys.stdin.isatty():
        raise DownloadError("Run --signed-links in an interactive terminal for hidden URL input.")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    opener = urllib.request.build_opener(NoRedirect())
    for record in pending:
        print("Downloading " + record["filename"])
        if args.signed_links:
            url = getpass.getpass("Temporary download URL for this archive (not saved): ").strip()
        else:
            url = obtain_url(opener, manifest["api_endpoint"], record, authorization)
        download(opener, url, record, args.out_dir / record["filename"])
        print(record["filename"] + ": SHA256 OK")
    receipt = {
        "download_verification_passed": True,
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "release": manifest["release"], "genome": manifest["genome"],
        "historical_release_confirmed": False,
        "selection_status": manifest["release_selection"]["status"],
        "source_manifest_sha256": sha256(MANIFEST),
        "files": [{"filename": r["filename"], "sha256": r["sha256"]} for r in manifest["files"]],
    }
    receipt_path = args.out_dir / ("download_verification_" + uuid.uuid4().hex + ".json")
    with receipt_path.open("x") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    print("Verification receipt: " + str(receipt_path))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (DownloadError, ValueError) as error:
        print("Download failed: " + str(error), file=sys.stderr)
        sys.exit(1)
    except (EOFError, KeyboardInterrupt):
        print("Download cancelled.", file=sys.stderr)
        sys.exit(1)

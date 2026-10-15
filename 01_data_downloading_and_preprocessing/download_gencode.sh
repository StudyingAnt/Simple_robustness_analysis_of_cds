#!/usr/bin/env bash
# Download the specified GENCODE protein-coding transcript FASTA and unpack it.
set -euo pipefail

release=40
out_dir=''
while (($#)); do
    case "$1" in
        --release) release="${2:?--release needs a value}"; shift 2 ;;
        --out-dir) out_dir="${2:?--out-dir needs a directory}"; shift 2 ;;
        --help|-h)
            echo 'Usage: bash download_gencode.sh --out-dir DIRECTORY [--release 40]'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
if [[ -z "$out_dir" || ! "$release" =~ ^[1-9][0-9]*$ ]]; then
    echo 'Specify --out-dir DIRECTORY and a positive integer release.' >&2
    exit 2
fi

mkdir -p -- "$out_dir"
out_dir=$(cd -- "$out_dir" && pwd)
filename="gencode.v${release}.pc_transcripts.fa"
compressed="$out_dir/$filename.gz"
fasta="$out_dir/$filename"
manifest="$out_dir/$filename.download.json"
for file in "$compressed" "$fasta" "$manifest"; do
    if [[ -e "$file" ]]; then
        echo "Refusing to overwrite: $file" >&2
        exit 1
    fi
done

url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${release}/${filename}.gz"
started_at=$(date -u +%Y-%m-%dT%H:%M:%SZ)
staging=$(mktemp -d "$out_dir/.gencode-download.XXXXXX")
trap 'rm -rf -- "$staging"' EXIT
echo "Downloading GENCODE release $release"
curl --fail --location --silent --show-error --retry 3 --connect-timeout 20 --max-time 600 \
    --output "$staging/$filename.gz" "$url"
gzip -t -- "$staging/$filename.gz"
gzip -dc -- "$staging/$filename.gz" > "$staging/$filename"

python3 - "$staging" "$filename" "$release" "$url" "$started_at" "$0" <<'PY'
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

stage, name, release, url, started, script = sys.argv[1:]
script = Path(script).resolve()
def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

files = []
for filename in [name + '.gz', name]:
    path = Path(stage) / filename
    entry = {'filename': filename, 'bytes': path.stat().st_size, 'sha256': sha256(path)}
    if filename == name:
        with path.open('rb') as stream:
            entry['fasta_records'] = sum(line.startswith(b'>') for line in stream)
    files.append(entry)
commit = subprocess.check_output(['git', '-C', str(script.parent), 'rev-parse', 'HEAD'], text=True).strip()
record = {'source_url': url, 'gencode_release': int(release), 'started_at_utc': started,
          'completed_at_utc': datetime.now(timezone.utc).isoformat(),
          'gzip_integrity_passed': True, 'script_sha256': sha256(script),
          'analysis_commit': commit, 'files': files}
(Path(stage) / (name + '.download.json')).write_text(json.dumps(record, indent=2) + '\n')
PY

mv -- "$staging/$filename.gz" "$compressed"
mv -- "$staging/$filename" "$fasta"
mv -- "$staging/$filename.download.json" "$manifest"
echo "Downloaded and unpacked: $fasta"
echo "Provenance: $manifest"

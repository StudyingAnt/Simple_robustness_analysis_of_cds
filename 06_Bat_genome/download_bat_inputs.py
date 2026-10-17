"""Download the manuscript's Bat1K assemblies and annotation snapshots."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import time
import urllib.request


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def download(out_dir, sources_file, annotations_only=False):
    out_dir, sources_file = Path(out_dir).resolve(), Path(sources_file).resolve()
    if out_dir.exists():
        raise FileExistsError(f'Choose a new output directory: {out_dir}')
    sources = json.loads(sources_file.read_text())
    assets = [a for a in sources['assets'] if not annotations_only or a['role'] != 'genome']
    if not assets or len({a['filename'] for a in assets}) != len(assets):
        raise ValueError('Empty or duplicate source assets')
    for item in assets:
        if Path(item['filename']).name != item['filename'] or item['expected_bytes'] < 1:
            raise ValueError('Invalid source asset')
    out_dir.parent.mkdir(exist_ok=True, parents=True)
    results = []
    with tempfile.TemporaryDirectory(prefix='.bat1k-download-', dir=out_dir.parent) as temporary:
        stage = Path(temporary) / 'raw'
        stage.mkdir()
        for asset in assets:
            path = stage / asset['filename']
            digest, total = hashlib.sha256(), 0
            with urllib.request.urlopen(asset['url'], timeout=60) as response, path.open('xb') as stream:
                for block in iter(lambda: response.read(1024 * 1024), b''):
                    stream.write(block)
                    digest.update(block)
                    total += len(block)
            if total != asset['expected_bytes']:
                raise ValueError(f'Size mismatch: {path.name}')
            actual = digest.hexdigest()
            if asset.get('sha256') and actual != asset['sha256']:
                raise ValueError(f'Reference SHA256 mismatch: {path.name}')
            gzip_checked = False
            if path.suffix == '.gz':
                with gzip.open(path, 'rb') as stream:
                    for _ in iter(lambda: stream.read(1024 * 1024), b''):
                        pass
                gzip_checked = True
            results.append(dict(filename=path.name, url=asset['url'], bytes=total,
                                sha256=actual, reference_sha256_checked=bool(asset.get('sha256')),
                                gzip_integrity_checked=gzip_checked, role=asset['role']))
            print(f'Downloaded and checked {path.name}: {total} bytes', flush=True)
        report = dict(download_passed=True, annotations_only=annotations_only, outputs=results,
                      sources_sha256=sha256(sources_file), script_sha256=sha256(__file__),
                      completed_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                      note='Raw files only. No CDS extraction or RUMS simulation performed.')
        (stage / 'download_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
        if out_dir.exists():
            raise FileExistsError(f'Output appeared during download: {out_dir}')
        stage.rename(out_dir)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True, help='New directory for raw files')
    parser.add_argument('--sources', type=Path, default=Path(__file__).with_name('bat1k_pilot_sources.json'))
    parser.add_argument('--annotations-only', action='store_true', help='Omit the two large genome FASTAs')
    args = parser.parse_args(argv)
    report = download(args.out_dir, args.sources, args.annotations_only)
    print(json.dumps(dict(download_passed=report['download_passed'], assets=len(report['outputs']),
                         annotations_only=report['annotations_only']), indent=2))
    print('Manifest:', args.out_dir / 'download_manifest.json')


if __name__ == '__main__':
    main()

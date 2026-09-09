import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('runs', nargs='+')
parser.add_argument('--workers', type=int, choices=[1, 2], default=2)
args = parser.parse_args()
rows = list(csv.DictReader(open('config/heat_downloads.tsv'), delimiter='\t'))
if set(args.runs) - {r['run'] for r in rows}:
    raise ValueError('Run absent from the reviewed download list')
def download(row):
    path = Path(row['path'])
    expected = int(row['bytes'])
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.stat().st_size != expected:
            raise ValueError(f'Existing input has wrong size: {path}')
        return
    if shutil.disk_usage(path.parent).free < expected * 1.2:
        raise ValueError('Insufficient free disk including 20% headroom')
    temporary = path.with_suffix(path.suffix + '.part')
    subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error',
                    '--retry', '2', '--connect-timeout', '20', '--max-time', '3600',
                    '--speed-limit', '1000', '--speed-time', '120', '--continue-at', '-',
                    '--output', str(temporary), row['url']], check=True)
    if temporary.stat().st_size != expected:
        raise ValueError(f'Incomplete transfer: {temporary}')
    temporary.rename(path)
    print(f"{row['run']} mate {row['mate']}: {expected} bytes", flush=True)

with ThreadPoolExecutor(max_workers=args.workers) as executor:
    list(executor.map(download, [r for r in rows if r['run'] in args.runs]))

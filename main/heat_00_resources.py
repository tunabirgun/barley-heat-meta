import argparse
import csv
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='Filenames from config/resource_downloads.tsv')
args = parser.parse_args()
with open('config/resource_downloads.tsv', encoding='utf-8') as stream:
    rows = list(csv.DictReader(stream, delimiter='\t'))
catalog = {Path(row['path']).name: row for row in rows}
if len(catalog) != len(rows) or set(args.files) - catalog.keys():
    raise ValueError('Unknown or duplicate resource filename')
for name in args.files:
    row = catalog[name]
    path = Path(row['path'])
    if path.is_absolute() or '..' in path.parts:
        raise ValueError('Resource destination must stay inside the repository')
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.part')
    subprocess.run(['curl', '--fail', '--location', '--retry', '2', '--continue-at', '-',
                    '--output', str(temporary), row['url']], check=True)
    if not temporary.stat().st_size:
        raise ValueError('Empty download')
    temporary.rename(path)

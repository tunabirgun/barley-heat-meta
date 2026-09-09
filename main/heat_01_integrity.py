import argparse
import csv
import gzip
import hashlib
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('runs', nargs='+')
parser.add_argument('--download-only', action='store_true')
args = parser.parse_args()
downloads = list(csv.DictReader(open('config/heat_downloads.tsv'), delimiter='\t'))
output = Path('results/heat/download_verification' if args.download_only else 'results/heat/integrity')
output.mkdir(parents=True, exist_ok=True)
sequence_symbols = b'ACGTRYSWKMBDHVNacgtryswkmbdhvn'
quality_symbols = bytes(range(33, 127))
for run in args.runs:
    rows = sorted((r for r in downloads if r['run'] == run), key=lambda r: int(r['mate']))
    if len(rows) != 2 or [r['mate'] for r in rows] != ['1', '2']:
        raise ValueError(f'Expected two distinct mates: {run}')
    for row in rows:
        path = Path(row['path'])
        if path.stat().st_size != int(row['bytes']):
            raise ValueError(f'Byte mismatch: {path}')
        with path.open('rb') as stream:
            observed = hashlib.file_digest(stream, 'md5').hexdigest()
        if observed != row['md5'].lower():
            raise ValueError(f'MD5 mismatch: {path}')
    if args.download_only:
        with open(output / f'{run}.tsv', 'w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=[*rows[0], 'status'], delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(dict(row, status='size_and_md5_verified') for row in rows)
        print(f'{run}: size and MD5 verified; full FASTQ integrity not checked', flush=True)
        continue
    records = 0
    with gzip.open(rows[0]['path'], 'rb') as first, gzip.open(rows[1]['path'], 'rb') as second:
        while True:
            pair = [[stream.readline() for _ in range(4)] for stream in (first, second)]
            if not any(line for record in pair for line in record):
                break
            identifiers = []
            for record in pair:
                header, sequence, plus, quality = [line.rstrip(b'\r\n') for line in record]
                if not header.startswith(b'@') or not plus.startswith(b'+') or not sequence or len(sequence) != len(quality):
                    raise ValueError(f'Invalid four-line FASTQ record: {run}, record {records + 1}')
                if sequence.translate(None, sequence_symbols) or quality.translate(None, quality_symbols):
                    raise ValueError(f'Invalid sequence or quality encoding: {run}')
                token = header.split()[0]
                identifiers.append(token[:-2] if token.endswith((b'/1', b'/2')) else token)
            if identifiers[0] != identifiers[1]:
                raise ValueError(f'Mate identity mismatch: {run}, record {records + 1}')
            records += 1
    if records == 0:
        raise ValueError(f'Empty FASTQ: {run}')
    with open(output / f'{run}.tsv', 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=[*rows[0], 'read_pairs', 'status'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(dict(row, read_pairs=records, status='passed') for row in rows)
    print(f'{run}: integrity passed; {records} read pairs', flush=True)

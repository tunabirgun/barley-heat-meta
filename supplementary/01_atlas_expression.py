import csv
import gzip
import math
import re
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

samples = list(csv.DictReader(open(sys.argv[1] if len(sys.argv)>1 else 'config/atlas_panel_samples.tsv'), delimiter='\t'))
output = Path(sys.argv[2] if len(sys.argv)>2 else 'results/atlas_panel')
files = {row['quant_file']: row['run'] for row in samples}
if not samples or len(files)!=len(samples) or len({row['run'] for row in samples})!=len(samples):
    raise ValueError('Empty or duplicated atlas sample identifiers')
lengths = defaultdict(int)
with gzip.open('data/atlas/PanBaRT20_renamed.gff3.gz','rt') as annotation:
    for line in annotation:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields)==9 and fields[2]=='exon':
            attributes = dict(value.split('=',1) for value in fields[8].strip(';').split(';'))
            for transcript in attributes['Parent'].split(','):
                lengths[transcript] += int(fields[4])-int(fields[3])+1
if not lengths:
    raise ValueError('Atlas annotation contains no transcript exon lengths')
expression = {}
with tarfile.open('data/atlas/panbart20_tpm_quant_files.tar.gz', 'r|gz') as archive:
    for item in archive:
        name = Path(item.name).name
        if name not in files:
            continue
        if not item.isfile() or files[name] in expression:
            raise ValueError('Invalid or duplicated selected archive member')
        values = defaultdict(float)
        observed = set()
        rows = csv.DictReader((line.decode() for line in archive.extractfile(item)), delimiter='\t')
        if rows.fieldnames != ['Name','Length','EffectiveLength','TPM','NumReads']:
            raise ValueError('Unexpected Salmon quantification header')
        for row in rows:
            transcript = re.sub(r'^(chr(?:[1-7]H|Un))(\d+)(\.\d+)$',r'PanBaRT20_\1G\2\3',row['Name'])
            if transcript in observed or transcript not in lengths or int(row['Length'])!=lengths[transcript]:
                raise ValueError('Transcript identity or length differs from the atlas annotation')
            numbers = [float(row[key]) for key in ('TPM','EffectiveLength','NumReads')]
            if any(not math.isfinite(value) or value<0 for value in numbers):
                raise ValueError('Invalid Salmon quantification value')
            values[transcript.rsplit('.',1)[0]] += numbers[0]
            observed.add(transcript)
        if observed!=lengths.keys() or not math.isclose(math.fsum(values.values()),1e6,rel_tol=1e-5):
            raise ValueError('Transcript universe or TPM total differs from the expected library')
        expression[files[name]] = values
if set(expression)!={row['run'] for row in samples}:
    raise ValueError('Selected atlas libraries are missing')
runs = [row['run'] for row in samples]
genes = sorted(expression[runs[0]])
if not all(set(values)==set(genes) for values in expression.values()):
    raise ValueError('Atlas gene universes differ')
output.mkdir(parents=True,exist_ok=True)
with open(output/'gene_tpm.tsv','w',newline='') as stream:
    writer = csv.writer(stream,delimiter='\t',lineterminator='\n')
    writer.writerow(['pan_gene_id',*runs])
    writer.writerows([gene,*[expression[run][gene] for run in runs]] for gene in genes)

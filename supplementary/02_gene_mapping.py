import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

morex_targets, pan_targets, exact = defaultdict(set), defaultdict(set), []
for filename, morex_query in [('morex_to_pan.tmap.gz',True),('pan_to_morex.tmap.gz',False)]:
    matches = set()
    with gzip.open('data/atlas/'+filename,'rt') as source:
        for row in csv.DictReader(source,delimiter='\t'):
            if row['class_code'] not in {'=','c','k','m','n','j','e','o'} or row['ref_id']=='-':
                continue
            morex = re.match(r'HORVU\.MOREX\.r3\.(?:[1-7]H|Un)G\d+',row['qry_id'] if morex_query else row['ref_id'])[0]
            pan = (row['ref_id'] if morex_query else row['qry_id']).rsplit('.',1)[0]
            morex_targets[morex].add(pan)
            pan_targets[pan].add(morex)
            if row['class_code']=='=':
                matches.add((morex,pan))
    exact.append(matches)
reciprocal = exact[0] & exact[1]
Path('results/atlas').mkdir(parents=True,exist_ok=True)
with open('results/atlas/gene_mapping.tsv','w',newline='') as output:
    writer = csv.writer(output,delimiter='\t',lineterminator='\n')
    writer.writerow(['gene_id','pan_gene_id','display_mapping'])
    for gene in sorted(morex_targets):
        for pan in sorted(morex_targets[gene]):
            keep = (gene,pan) in reciprocal and len(morex_targets[gene])==len(pan_targets[pan])==1
            writer.writerow([gene,pan,'reciprocal_unique' if keep else 'ambiguous_or_nonexact'])

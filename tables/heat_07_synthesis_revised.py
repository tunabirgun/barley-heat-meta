import csv
import json
from pathlib import Path

import numpy as np

root = Path('.')
source = root / 'results/heat/synthesis_revised'
out = root / 'results/heat/presentation'
tables = out / 'tables'
tables.mkdir(parents=True, exist_ok=True)

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(path, rows, fields):
    if path.exists():
        raise FileExistsError(path)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t'); w.writeheader(); w.writerows(rows)

pc = read(source / 'partial_conjunction.tsv')
re = {x['gene_id']: x for x in read(source / 'random_effects.tsv')}
fam = read(source / 'family_effects.tsv')
families = list(dict.fromkeys(x['publication_family_id'] for x in read(root/'config/heat_models.tsv')))
by_family = {(x['gene_id'], x['family']): x for x in fam}
all_q = f'q_at_least_{len(families)}'
candidates = [x for x in pc if float(x[all_q]) <= .05]
candidates.sort(key=lambda x: (-abs(float(re[x['gene_id']]['reml_estimate'])) if re[x['gene_id']]['reml_estimate'] != 'NA' else float('inf'), x['gene_id']))
descriptions = {x['gene_id']: x for x in read(root/'config/gene_descriptions.tsv')}
heat_genes = {x['gene_id'] for x in read(root/'config/gene_go.tsv') if x['pathway_id'] in {'GO:0009408','GO:0034605'}}
mapping = read(root/'results/atlas/gene_mapping.tsv')
rows = []
for rank, c in enumerate(candidates, 1):
    g = c['gene_id']; matches = [x for x in mapping if x['gene_id'] == g]
    unique = [x for x in matches if x['display_mapping'] == 'reciprocal_unique']
    if len(unique) > 1:
        raise ValueError(f'Duplicate reciprocal mapping: {g}')
    status = 'mapped' if unique else 'ambiguous' if matches else 'unmapped'
    d = descriptions[g]
    row = dict(gene_id=g, rank=rank, direction=c['direction'], partial_conjunction_q_all=c[all_q],
        **{k:v for k,v in re[g].items() if k!='gene_id'}, heat_go_annotated=g in heat_genes,
        annotation=next((d[k] for k in ('description','gene_name') if d.get(k) not in (None,'','NA')), 'not annotated'),
        atlas_mapping=status, pan_gene_id=unique[0]['pan_gene_id'] if unique else 'NA')
    for family in families:
        f = by_family[g,family]
        row.update({f'{family}_log2FoldChange':f['estimate'], f'{family}_se':f['se'], f'{family}_q_within_family':f['q_within_family']})
    rows.append(row)
fields = list(rows[0]) if rows else ['gene_id','rank','direction','partial_conjunction_q_all']
write(tables/'candidates.tsv', rows, fields)
selected = {x['gene_id'] for x in rows}
effects = {}
for path in [root/'results/heat/inference/ordinary_results.tsv', root/'results/heat/shrink/stabilized_display.tsv']:
    with path.open(encoding='utf-8-sig', newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['gene_id'] in selected:
                effects.setdefault((r['gene_id'],r['contrast_id']), {}).update(r)
values = list(effects.values())
if values:
    write(tables/'candidate_contrast_values.tsv', values, list(values[0]))
samples = read(root/'config/atlas_panel_samples.tsv')
strains = list(dict.fromkeys(x['strain'] for x in samples))
assert len(strains)==20 and all(sum(x['strain']==s for x in samples)==3 for s in strains)
pan_to_genes = {}
for r in rows:
    if r['atlas_mapping']=='mapped':pan_to_genes.setdefault(r['pan_gene_id'], []).append(r['gene_id'])
atlas = []
with (root/'results/atlas_panel/gene_tpm.tsv').open(encoding='utf-8-sig', newline='') as f:
    reader=csv.reader(f, delimiter='\t'); header=next(reader)
    for r in reader:
        if r[0] not in pan_to_genes:continue
        tpm={run:float(v) for run,v in zip(header[1:],r[1:])}
        assert all(np.isfinite(v) and v>=0 for v in tpm.values())
        for g in pan_to_genes[r[0]]:
            for s in samples:
                atlas.append(dict(gene_id=g,pan_gene_id=r[0],strain=s['strain'],run=s['run'],tpm=tpm[s['run']],log2_tpm_plus1=float(np.log2(tpm[s['run']]+1))))
write(tables/'atlas_replicates.tsv', atlas, ['gene_id','pan_gene_id','strain','run','tpm','log2_tpm_plus1'])
summary=[]
for g in [x['gene_id'] for x in rows]:
    for s in strains:
        v=[x['log2_tpm_plus1'] for x in atlas if x['gene_id']==g and x['strain']==s]
        summary.append(dict(gene_id=g,strain=s,n=len(v),mean_log2_tpm_plus1=float(np.mean(v)) if v else 'NA',sd_log2_tpm_plus1=float(np.std(v,ddof=1)) if v else 'NA'))
write(tables/'atlas_summary.tsv', summary, ['gene_id','strain','n','mean_log2_tpm_plus1','sd_log2_tpm_plus1'])
audit=dict(candidates=len(rows),up=sum(x['direction']=='up' for x in rows),down=sum(x['direction']=='down' for x in rows),
    mapped=sum(x.get('atlas_mapping')=='mapped' for x in rows),heat_annotated=sum(x.get('heat_go_annotated',False) for x in rows),
    median_I2=float(np.median([float(x['reml_I2']) for x in rows if x['reml_I2']!='NA'])) if rows else None,
    pooled_q05=sum(x.get('reml_q')!='NA' and float(x['reml_q'])<=.05 for x in rows),
    covariance='upper_bound_not_estimated',tier='conditional_exploratory')
(out/'summary.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit))

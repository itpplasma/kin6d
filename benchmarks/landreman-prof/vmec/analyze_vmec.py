"""Medians of the VMEC++ ladder job: python3 analyze_vmec.py OUTDIR  (OUTDIR contains vmec/ and points.txt)."""
import json, re, sys
from collections import defaultdict
from pathlib import Path
from statistics import median
out = Path(sys.argv[1]); keys = ['B_sample_rms_relative=', 'p_sample_rms_relative=', 'surface_label_max_relative=', 'force_sample_rms_relative=']
spec = {l.split()[0]: l.split()[1:] for l in (out / 'points.txt').read_text().splitlines() if l.strip()}
rows = defaultdict(lambda: defaultdict(list))
for d in sorted((out / 'vmec').glob('*-r[0-9]')):
    tag, t, r = d.name.rsplit('-', 2)
    try: m = json.loads((d / 'metrics.json').read_text())
    except Exception: m = {}
    c = (d / 'compare.log').read_text() if (d / 'compare.log').exists() else ''
    v = rows[(tag, int(t[1:]))]; v['n'].append(1)
    v['native'].append(m.get('native_solve_wall_s')); v['wf'].append(m.get('solve_wall_s'))
    for k, x in zip(['fsqr', 'fsqz', 'fsql'], [m.get(k) for k in ('fsqr', 'fsqz', 'fsql')]): v[k].append(x)
    for k in keys:
        mm = re.search(re.escape(k) + r'\s*(\S+)', c); v[k].append(float(mm.group(1)) if mm else None)
def med(x):
    x = [a for a in x if a is not None]; return median(x) if x else None
def f(x, fmt='{:.4g}'): return '-' if x is None else fmt.format(x)
print('|pt|mpol|NS_ARRAY|FTOL_ARRAY|t|ok/3|native s|workflow s|B %|p %|label %|max %|force rel|fsqr|\n|---|---|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|')
res = []
for (tag, t), v in sorted(rows.items()):
    e = [med(v[k]) for k in keys]; ok = sum(a is not None for a in v['native'])
    mx = max(e[:3]) * 100 if None not in e[:3] else None
    s = spec.get(tag, ['?'] * 4)
    print(f"|{tag}|{s[0]}|{s[1]}|{s[2]}|{t}|{ok}/{len(v['n'])}|{f(med(v['native']))}|{f(med(v['wf']))}|" +
          '|'.join(f(x * 100 if x is not None else None, '{:.3f}') for x in e[:3]) + f"|{f(mx, '{:.3f}')}|{f(e[3])}|{f(med(v['fsqr']))}|")
    if mx is not None and med(v['native']): res.append((t, med(v['native']), tag, mx))
print()
for target in (1.0, 0.5, 0.25):
    for t in (1, 8, 16):
        c = sorted(x for x in res if x[0] == t and x[3] <= target)
        print(f'- target {target}% t={t}: ' + ('; '.join(f'{n} {w:.4g}s max {e:.3f}%' for _, w, n, e in c[:3]) or 'none'))
fl = out / 'failures.txt'; print('\nfailures:', fl.read_text().strip() if fl.exists() else 'none')

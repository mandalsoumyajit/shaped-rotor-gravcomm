"""Archive summary-only markers before regenerating omitted FEA solver fields."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
for name in ['three_link_designs_2026_09_20','tungsten_comparison_2026_09_20']:
 for p in (ROOT/'fea/results'/name).rglob('summary.json'):
  if all((p.parent/n).exists() for n in ['model.inp','model.dat','model.frd','fields.npz']):continue
  backup=p.with_name('summary_published.json')
  if backup.exists():raise RuntimeError(f'Backup already exists: {backup}')
  p.rename(backup);print('Archived summary marker:',p.relative_to(ROOT))

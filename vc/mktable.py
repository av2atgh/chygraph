"""results/compare.json + results/synthetic.json -> results/compare_table.tex"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parent / 'results'
out = []

rows = json.loads((RES / 'compare.json').read_text()) if (RES / 'compare.json').exists() else []
if rows:
    out.append(r'\paragraph{Degree-correlated ensemble} ($\tau=2.5$, seed 0; exact from HiGHS, '
               r'time-limited; ``node'' is the node-level min-sum, ``complex'' keeps SBM blocks of '
               r'$\le128$ vertices as exactly solved regions; $E$ is the Bethe energy, $|C|$ the '
               r'decimated cover).')
    out.append(r'\begin{center}\small\begin{tabular}{rrrrrrrrrrr}\toprule')
    out.append(r'$r$ & $n$ & exact & leaf rem. & core & $E_{\rm node}$ & $|C|_{\rm node}$ & conv. & '
               r'$E_{\rm complex}$ & $|C|_{\rm complex}$ & conv.\\ \midrule')
    for x in sorted(rows, key=lambda x: (x['n'], x['r'])):
        bp = {b['smax']: b for b in x['bp']}
        a, b = bp[0], bp[max(bp)]
        ex = f"{x['exact']:.0f}" if x['exact'] is not None else '--'
        if x['exact'] is not None and x['exact_gap'] and x['exact_gap'] > 1e-6:
            ex += f"$^*$"
        conv = lambda z: 'yes' if z['resid'] < 1e-6 else 'no'
        out.append(f"{x['r']:.1f} & {x['n']} & {ex} & {x['leaf']} & {x['core']} & "
                   f"{a['bethe']:.1f} & {a['cover']} & {conv(a)} & {b['bethe']:.1f} & {b['cover']} & {conv(b)}\\\\")
    out.append(r'\bottomrule\end{tabular}\end{center}')
    out.append(r'$^*$: MILP stopped at the time limit; value is the incumbent, not proven optimal.')

syn = sum((json.loads(p.read_text()) for p in sorted(RES.glob('synthetic*.json'))), [])
if syn:
    out.append(r'\paragraph{Planted dense blocks} (150 blocks of 5--40 vertices, internal density '
               r'$\rho$, plus $c$ random inter-block links per vertex).')
    out.append(r'\begin{center}\small\begin{tabular}{rrrrrrrrrr}\toprule')
    out.append(r'$\rho$ & $c$ & $n$ & exact & core & regions & $E_{\rm node}$ & $|C|_{\rm node}$ & '
               r'$E_{\rm complex}$ & $|C|_{\rm complex}$\\ \midrule')
    keys = sorted({(s['rho'], s['c']) for s in syn}, reverse=True)
    for rho, c in keys:
        sel = [s for s in syn if s['rho'] == rho and s['c'] == c]
        node = [s for s in sel if s['smax'] == 0][0]
        for s in [s for s in sel if s['smax'] > 0]:
            out.append(f"{rho} & {c} & {s['n']} & {s['exact']:.0f} & {s['core']} & {s['label']} ({s['regions']}) & "
                       f"{node['bethe']:.1f} & {node['cover']} & {s['bethe']:.1f} & {s['cover']}\\\\")
    out.append(r'\bottomrule\end{tabular}\end{center}')

(RES / 'compare_table.tex').write_text('\n'.join(out) + '\n')
print('\n'.join(out))

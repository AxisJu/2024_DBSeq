def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    from pathlib import Path
    import os
    import json
    import numpy as np
    import pandas as pd
    from scipy import stats
    PACKAGE = Path(os.environ['DBSEQ_OUTPUT_ROOT'])
    code_base = Path(__file__).resolve().parents[2]
    if PACKAGE.resolve() == code_base or code_base in PACKAGE.resolve().parents:
        raise ValueError('Output must be outside the code package')
    DATA = Path(os.environ['DBSEQ_DATA_ROOT'])
    WORK = PACKAGE
    FIG = PACKAGE / 'figures/extended data figure 13'
    SRC = PACKAGE / 'sourcedata/extended data figure 13'
    AUDIT = WORK / 'cache/extended_data'
    for p in (FIG, SRC, AUDIT):
        p.mkdir(parents=True, exist_ok=True)
    METRICS = ['FC', 'ED', 'SC', 'SC2', 'SPE', 'DE', 'NE', 'CMY']
    BINS = ['2-4', '5-7', '8-10', '11-13', '14-16', '17-19', '20+']
    DIV = {'CNU': '#9bd5f4', 'Isocortex': '#18a799', 'TH': '#e0abce', 'HY': '#e74538', 'HPF': '#84c551'}

    def frontier(frame, metric, quantile):
        use = frame[[metric, 'expr']].dropna().sort_values(metric)
        if use.empty:
            return np.nan
        ix = np.searchsorted(np.cumsum(use.expr) / use.expr.sum(), quantile)
        return float(use.iloc[ix][metric])

    def scatter(ax, frame, metric, partial):
        x = frame.commonness.to_numpy()
        y = frame[metric].to_numpy()
        if partial:
            rx, ry, rc = map(stats.rankdata, (x, y, frame.ED.to_numpy()))
            rx -= np.polyval(np.polyfit(rc, rx, 1), rc)
            ry -= np.polyval(np.polyfit(rc, ry, 1), rc)
            rho = stats.pearsonr(rx, ry).statistic
            p = 2 * stats.t.sf(abs(rho) * np.sqrt((len(x) - 3) / (1 - rho * rho)), len(x) - 3)
        else:
            rho, p = stats.spearmanr(x, y)
        return {'rho': float(rho), 'p': float(p), 'n': len(x)}

    def main():
        membership = inputs['membership_1'][str(AUDIT / 'ed13_membership.csv')]
        weights = inputs['weights_2'][str(AUDIT / 'ed13_frontier_membership.csv')]
        imaging = inputs['imaging_3'][str(DATA / 'results/Table/imaging_metrics_8metrics_12_regions.csv')]
        mapping = {'CP': 'STRd', 'LSc': 'LSX', 'TRS': 'LSX', 'CM': 'ILM', 'PF': 'ILM', 'RE': 'MTN', 'GPe': 'PALd', 'MH': 'EPI', 'LH': 'LZ', 'RT': 'RT', 'RSPd': 'RSPd', 'DG': 'DG'}
        fc = inputs['fc_4'][str(DATA / 'data/derivatives/dbseq/mouse_BOLD_fc.csv')]
        imaging['FC'] = [float(fc.loc['ATN', mapping[r]]) for r in imaging.Region]
        imaging.loc[imaging.Region.eq('MH'), 'Division'] = 'TH'
        weighted = weights.merge(imaging, on='Region', validate='many_to_one')
        rows = []
        for (gene, n), group in weighted.groupby(['genes', 'n_neurons']):
            row = {'gene': gene, 'commonness': n}
            row.update({m: frontier(group, m, 0.9 if m == 'ED' else 0.1) for m in METRICS})
            rows.append(row)
        gene = pd.DataFrame(rows)
        gene.assign(unit='neuron subtype count; ED micrometre; weighted 10th metric or 90th ED percentile').to_csv(SRC / 'ed13b.csv', index=False)
        b_stats = {}
        for m in METRICS[:4]:
            x = gene.commonness
            y = gene[m]
            rho, p = stats.spearmanr(x, y)
            b_stats[m] = {'rho': float(rho), 'p': float(p)}
        membership['bin'] = pd.cut(membership.n_neurons, [1, 4, 7, 10, 13, 16, 19, np.inf], labels=BINS)
        region = membership.groupby('Region').n_neurons.mean().rename('commonness').reset_index().merge(imaging, on='Region', validate='one_to_one')
        counts = membership.groupby(['bin', 'Region'], observed=True).size().rename('count').reset_index()
        counts['proportion'] = counts['count'] / counts.groupby('bin', observed=True)['count'].transform('sum')
        membership.merge(counts, on=['bin', 'Region']).assign(unit='gene-neuronal subtype occurrences; proportion fraction').to_csv(SRC / 'ed13c.csv', index=False)
        order = region.sort_values('commonness').Region.tolist()
        matrix = counts.pivot(index='bin', columns='Region', values='proportion').reindex(index=BINS, columns=order).fillna(0)
        report = {'b': b_stats, 'd': {}, 'e': {}, 'definition': 'Regional commonness weighted by gene-celltype occurrences, not number of individual nuclei. Gene frontier uses positive CosSim deviation weights and 10th/90th percentile as original.'}
        for panel, metrics, partial in [('d', METRICS[:4], False), ('e', METRICS[4:], True)]:
            cols = ['Region', 'commonness', 'Division'] + list(dict.fromkeys(metrics + ['ED']))
            region[cols].assign(unit='neuronal subtype mean count; ED micrometre; network metric native scale').to_csv(SRC / f'ed13{panel}.csv', index=False)
            report[panel] = {m: scatter(None, region, m, partial) for m in metrics}
        report['differences'] = 'Run ed13bcde.R after this entry to retain original quadratic LOESS display. Panel d FC uses original LH-to-LZ mapping and matches final rho=0.01. Panel e NE partial rho=-0.739 differs from final -0.59 but its P=0.009312 matches; raw NE matrix agrees. No observations manually excluded.'
        (AUDIT / 'ed13-audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(f'Extended Data Figure 13 sources: {len(gene)} genes, {len(region)} regions; run ed13bcde.R')
    main()
    return locals()

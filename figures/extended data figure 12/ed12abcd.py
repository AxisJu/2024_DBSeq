def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    from pathlib import Path
    import os
    import json
    import numpy as np
    import pandas as pd
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    PACKAGE = Path(os.environ['DBSEQ_OUTPUT_ROOT'])
    code_base = Path(__file__).resolve().parents[2]
    if PACKAGE.resolve() == code_base or code_base in PACKAGE.resolve().parents:
        raise ValueError('Output must be outside the code package')
    DATA = Path(os.environ['DBSEQ_DATA_ROOT'])
    WORK = PACKAGE
    FIG = PACKAGE / 'figures/extended data figure 12'
    SRC = PACKAGE / 'sourcedata/extended data figure 12'
    AUDIT = WORK / 'cache/extended_data'
    for p in (FIG, SRC, AUDIT):
        p.mkdir(parents=True, exist_ok=True)
    DIV = {'CNU': '#9bd5f4', 'Isocortex': '#18a799', 'TH': '#e0abce', 'HY': '#e74538', 'HPF': '#84c551'}
    METRICS = ['ED', 'FC', 'SC', 'SC2', 'SPE', 'DE', 'NE', 'CMY']

    def source(frame, panel):
        frame.to_csv(SRC / f'ed12{panel}.csv', index=False)

    def clean(term):
        return term.split(' (GO:')[0]

    def inputs():
        membership = inputs['membership_1'][str(AUDIT / 'ed12_regional_genes.csv')]
        regional = {r: set(g.Mouse_Gene.str.upper()) for r, g in membership.groupby('Region')}
        background = set(inputs['background_2'][str(AUDIT / 'ed12_background_genes.txt')].upper().split())
        gene_sets = inputs['gene_sets_3'][str(AUDIT / 'ed12_all_gene_sets.csv')]
        go = {row.pathway: set(row.genes.upper().split(';')) for row in gene_sets.itertuples()}
        imaging = inputs['imaging_4'][str(DATA / 'results/260813 DOLPHIN/Merged_Splicing_and_8Imaging_12_Regions.csv')]
        observations = []
        summaries = []
        atn = regional['ATN']
        for term, genes in go.items():
            genes &= background
            k = len(atn & genes)
            p = stats.hypergeom.sf(k - 1, len(background), len(background & genes), len(atn & background))
            frame = imaging[['Region', 'Division'] + METRICS].copy()
            frame['pathway'] = term
            frame.loc[frame.Region.eq('MH'), 'Division'] = 'TH'
            frame['TRG_hits'] = [len(regional[r] & genes) for r in frame.Region]
            frame['TRG_genes'] = [';'.join(sorted(regional[r] & genes)) for r in frame.Region]
            frame['ATN_TRG_hits'] = k
            frame['ATN_TRG_genes'] = ';'.join(sorted(atn & genes))
            frame['pathway_background_gene_count'] = len(background & genes)
            frame['background_gene_count'] = len(background)
            frame['ATN_TRG_total'] = len(atn)
            corr = stats.spearmanr(frame.TRG_hits, frame.SC2) if frame.TRG_hits.nunique() > 1 else (0.0, 1.0)
            summaries.append({'pathway': term, 'p_atn': float(p), 'rho': float(corr[0]), 'p_sc2': float(corr[1]), 'ATN_hits': k})
            observations.append(frame)
        summary = pd.DataFrame(summaries)
        summary['q_atn'] = multipletests(summary.p_atn, method='fdr_bh')[1]
        summary['q_sc2'] = 1.0
        eligible = summary.ATN_hits.ge(6)
        if eligible.any():
            summary.loc[eligible, 'q_sc2'] = multipletests(summary.loc[eligible, 'p_sc2'], method='fdr_bh')[1]
        return (pd.concat(observations, ignore_index=True), summary, go)

    def main():
        obs, summary, go = inputs()
        summary = summary.sort_values(['rho', 'pathway'], ascending=[False, True])
        obs['unit'] = 'TRG gene count; ED micrometre; connectivity metric native scale'
        source(obs[obs.pathway.isin(summary.loc[summary.ATN_hits.ge(6), 'pathway'])], 'b')
        eligible = summary[summary.ATN_hits.ge(6)]
        counts = [int(((eligible.rho > 0) & (eligible.q_sc2 < 0.05)).sum()), int(((eligible.rho > 0) & (eligible.q_sc2 >= 0.05)).sum()), int((eligible.rho < 0).sum())]
        fdr = summary[(summary.q_atn < 0.1) & summary.ATN_hits.ge(6)].copy()
        selected_obs = obs[obs.pathway.isin(fdr.pathway)]
        source(selected_obs[['pathway', 'Region', 'TRG_hits', 'TRG_genes', 'SC2', 'ATN_TRG_hits', 'pathway_background_gene_count', 'background_gene_count', 'ATN_TRG_total', 'unit']], 'a')
        retained = ['Aerobic Electron Transport Chain (GO:0019646)', 'Spliceosomal Complex Assembly (GO:0000245)', 'protein-DNA Complex Organization (GO:0071824)', 'Endoplasmic Reticulum To Golgi Vesicle-Mediated Transport (GO:0006888)', 'Regulation Of mRNA Splicing, Via Spliceosome (GO:0048024)', 'Regulation Of mRNA Stability (GO:0043488)', 'Protein Insertion Into Mitochondrial Membrane (GO:0051204)', 'Microtubule Depolymerization (GO:0007019)']
        assert set(retained).issubset(set(go))
        source(obs[obs.pathway.isin(retained)][['pathway', 'Region', 'Division', 'TRG_hits', 'TRG_genes', 'SC2', 'unit']], 'c')
        source(selected_obs[['pathway', 'Region', 'TRG_hits', 'TRG_genes'] + METRICS + ['unit']], 'd')
        statistics_dir = PACKAGE / 'statistics/extended data figure 12'
        statistics_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(statistics_dir / 'ed12_pathway_tests.csv', index=False)
        heat = []
        for term, frame in obs.groupby('pathway'):
            if term not in set(summary.loc[summary.ATN_hits.ge(6), 'pathway']):
                continue
            for metric in METRICS:
                r, p = stats.spearmanr(frame.TRG_hits, frame[metric]) if frame.TRG_hits.nunique() > 1 and frame[metric].nunique() > 1 else (0.0, 1.0)
                heat.append({'pathway': term, 'Metric': metric, 'rho': r, 'p': p})
        heat = pd.DataFrame(heat)
        heat['FDR'] = multipletests(heat.p, method='fdr_bh')[1]
        heat.to_csv(statistics_dir / 'ed12_all_metric_tests.csv', index=False)
        report = {'pathways': int(summary.ATN_hits.ge(6).sum()), 'full_ORA_family': len(summary), 'pie_counts': counts, 'ATN_FDR_under_0_1': len(fdr), 'scatter_selection': retained, 'selection_rule': 'Eight frozen formal display pathways; all regions retained; Displayed examples fixed independently of current significance; all regions retained', 'difference': 'Formal eight display examples restored; no observations removed for significance.', 'input_provenance': 'prepare_ed12.R regenerates membership and background from original QS files; ed12_all_gene_sets.csv contains the full frozen GO 2023 library, BH precedes display hit filtering; full-precision imaging values from splicing imaging table'}
        (AUDIT / 'ed12-audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(f'Extended Data Figure 12 sources (run ed12render.R): {len(summary)} pathways, {len(fdr)} FDR pathways')
    main()
    return locals()

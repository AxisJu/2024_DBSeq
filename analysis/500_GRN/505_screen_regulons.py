def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def bh(values):
        import numpy as np
        p = np.asarray(values, dtype=float)
        result = np.full(p.shape, np.nan)
        valid = np.flatnonzero(np.isfinite(p))
        order = valid[np.argsort(p[valid])]
        if len(order):
            result[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1])
        return result

    def donor_tests(auc, metadata, min_cells):
        import pandas as pd
        from scipy.stats import mannwhitneyu
        subsets = {'All': metadata.index, 'ExN': metadata.index[metadata.celltype_level1 == 'ExN'], 'InN': metadata.index[metadata.celltype_level1 == 'InN'], 'ATN': metadata.index[metadata.final_region == 'ATN']}
        for region in sorted(metadata.final_region.dropna().unique()):
            subsets[str(region)] = metadata.index[metadata.final_region == region]
        rows, observations = ([], [])
        for subset, ids in subsets.items():
            frame = auc.loc[ids].join(metadata.loc[ids, ['donor_id', 'Group_L2']])
            count = frame.groupby(['donor_id', 'Group_L2']).size()
            valid = count[count >= min_cells].index
            means = frame.groupby(['donor_id', 'Group_L2'])[auc.columns].mean().loc[valid]
            long = means.reset_index().melt(['donor_id', 'Group_L2'], var_name='regulon', value_name='AUC')
            long['subset'] = subset
            observations.append(long)
            for b, a in [('DBS_I', 'Sham_I'), ('Sham_I', 'Saline_I')]:
                for regulon in auc.columns:
                    x = means.loc[means.index.get_level_values('Group_L2') == b, regulon].dropna()
                    y = means.loc[means.index.get_level_values('Group_L2') == a, regulon].dropna()
                    p = mannwhitneyu(x, y, alternative='two-sided', method='asymptotic').pvalue if min(len(x), len(y)) >= 2 else float('nan')
                    rows.append(dict(subset=subset, regulon=regulon, comparison=b + '_vs_' + a, n_donors_b=len(x), n_donors_a=len(y), mean_b=x.mean(), mean_a=y.mean(), mean_difference=x.mean() - y.mean(), p_value=p))
        result = pd.DataFrame(rows)
        result['q_all_regulons_subsets_contrasts'] = bh(result.p_value)
        result['legacy_exploratory_p_lt_005'] = result.p_value < 0.05
        return (result, pd.concat(observations, ignore_index=True))

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('500_GRN')
        import pandas as pd
        from scipy.stats import fisher_exact, spearmanr
        resolve = lambda value, name: input_path(data, value) if value else out / name
        auc = inputs['auc_4'][str(resolve(a.auc, 'regulon_auc.csv'))]
        metadata = inputs['metadata_5'][str(resolve(a.metadata, 'regulon_metadata.csv'))]
        if auc.index.duplicated().any() or metadata.index.duplicated().any():
            raise ValueError('Duplicate cell IDs')
        metadata = metadata.loc[auc.index]
        if not {'donor_id', 'Group_L2', 'celltype_level1', 'final_region'} <= set(metadata.columns):
            raise ValueError('Missing donor or annotation columns')
        valid = metadata.donor_id.notna() & ~metadata.donor_id.isin(['unassigned', 'doublet', ''])
        metadata, auc = (metadata.loc[valid], auc.loc[valid])
        if (metadata.groupby('donor_id').Group_L2.nunique() > 1).any():
            raise ValueError('A donor occurs in multiple conditions; specify a paired design separately')
        tests, means = donor_tests(auc, metadata, a.min_cells)
        tests.to_csv(out / 'donor_regulon_tests.csv', index=False)
        means.to_csv(out / 'donor_regulon_auc.csv', index=False)
        wide = tests.pivot(index=['subset', 'regulon'], columns='comparison', values=['mean_difference', 'q_all_regulons_subsets_contrasts'])
        if all((name in wide['mean_difference'] for name in ['DBS_I_vs_Sham_I', 'Sham_I_vs_Saline_I'])):
            rescue = wide['mean_difference']['DBS_I_vs_Sham_I'] * wide['mean_difference']['Sham_I_vs_Saline_I'] < 0
            rescue &= wide['q_all_regulons_subsets_contrasts']['DBS_I_vs_Sham_I'] < 0.05
            rescue &= wide['q_all_regulons_subsets_contrasts']['Sham_I_vs_Saline_I'] < 0.05
            rescue.rename('reversed_at_BH_005').to_csv(out / 'exploratory_regulon_reversal.csv')
        if a.trgs:
            trgs = inputs['trgs_6'][str(input_path(data, a.trgs))]
            if trgs.gene.duplicated().any():
                raise ValueError('TRG background must have one row per tested gene')
            universe = set(trgs.gene)
            selected = set(trgs.loc[trgs.is_TRG.astype(str).str.lower().isin(['true', '1']), 'gene'])
            targets = inputs['targets_7'][str(resolve(a.targets, 'regulon_targets.csv'))]
            rows = []
            for tf, frame in targets.groupby('TF'):
                target = set(frame.target_gene) & universe
                hit = target & selected
                table = [[len(hit), len(target - selected)], [len(selected - target), len(universe - target - selected)]]
                rows.append(dict(TF=tf, n_targets_tested=len(target), n_TRGs=len(hit), n_background=len(universe), p_value=fisher_exact(table, alternative='greater').pvalue, target_TRGs=';'.join(sorted(hit))))
            enrich = pd.DataFrame(rows)
            enrich['q_all_TFs'] = bh(enrich.p_value)
            enrich.to_csv(out / 'regulon_TRG_overlap.csv', index=False)
        if a.imaging:
            imaging = inputs['imaging_8'][str(input_path(data, a.imaging))]
            if imaging.duplicated(['region', 'metric']).any():
                raise ValueError('Expected one imaging value per region and metric')
            regional = auc.join(metadata[['donor_id', 'Group_L2', 'final_region']]).groupby(['final_region', 'Group_L2', 'donor_id'])[auc.columns].mean().groupby(['final_region', 'Group_L2']).mean()
            rows = []
            for regulon in auc.columns:
                effects = regional[regulon].unstack('Group_L2')
                if not {'DBS_I', 'Sham_I'} <= set(effects.columns):
                    continue
                effects = (effects.DBS_I - effects.Sham_I).rename('auc_difference')
                for metric, frame in imaging.groupby('metric'):
                    joined = frame.set_index('region').join(effects).dropna(subset=['value', 'auc_difference'])
                    if len(joined) < 4:
                        continue
                    rho, pv = spearmanr(joined.value, joined.auc_difference)
                    rows.append(dict(regulon=regulon, metric=metric, n_regions=len(joined), rho=rho, p_value=pv))
            result = pd.DataFrame(rows)
            if len(result):
                result['q_all_regulons_metrics'] = bh(result.p_value)
            result.to_csv(out / 'regional_imaging_associations.csv', index=False)
    main()
    return locals()

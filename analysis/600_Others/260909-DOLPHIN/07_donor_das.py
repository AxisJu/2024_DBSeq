def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def load_events(se, mxe):
        import pandas as pd
        rows = []
        for kind, path, inclusion, exclusion in [('SE', se, ['junction12', 'junction23'], ['junction13']), ('MXE', mxe, ['junction12', 'junction24'], ['junction13', 'junction34'])]:
            table = inputs['table_3'][str(path)]
            for event, row in table.iterrows():
                inc, exc = ([str(row[k]) for k in inclusion], [str(row[k]) for k in exclusion])
                if any((k == 'nan' for k in inc + exc)):
                    raise ValueError('Incomplete event: ' + str(event))
                rows.append(dict(event_id=event, event_type=kind, gene_name=row.get('isoform1_gene_name', ''), inclusion=inc, exclusion=exc))
        return rows

    def bh(values):
        import numpy as np
        p = np.asarray(values, dtype=float)
        q = np.full(p.shape, np.nan)
        valid = np.flatnonzero(np.isfinite(p))
        order = valid[np.argsort(p[valid])]
        if len(order):
            q[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1])
        return q

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_4']
        data, out = roots('600_Others/260909-DOLPHIN')
        import numpy as np
        import pandas as pd
        from scipy.stats import permutation_test
        meta = inputs['meta_5'][str(input_path(data, a.metadata))]
        required = ['CB', 'donor_id', 'Group_L2', 'final_region']
        if not set(required) <= set(meta) or meta.CB.duplicated().any():
            raise ValueError('Unique cell metadata with biological donor IDs is required')
        if meta[required].isna().any().any() or meta.donor_id.isin(['unassigned', 'doublet', '']).any():
            raise ValueError('Unresolved donor or annotation')
        if (meta.groupby('donor_id').Group_L2.nunique() > 1).any():
            raise ValueError('Paired donors require a paired design; this entry uses independent donors')
        meta = meta[meta.Group_L2.isin(['DBS_I', 'Sham_I'])]
        events = inputs['events_6'][str(input_path(data, a.se_events))]
        junctions = {j for e in events for j in e['inclusion'] + e['exclusion']}
        pieces = []
        keys = ['final_region', 'Group_L2', 'donor_id', 'junction_id']
        for chunk in inputs['data_7'][str(input_path(data, a.reads))]:
            if chunk.sample_id.astype(str).str.contains('.aggr', regex=False).any():
                raise ValueError('Aggregated neighbor reads detected; supply original unaggregated reads')
            if 'read_origin' not in chunk or not chunk.read_origin.eq('original_STAR_unique').all():
                raise ValueError('Require read_origin=original_STAR_unique from 06_raw_junctions.py')
            if (chunk.reads < 0).any() or not np.equal(chunk.reads, np.round(chunk.reads)).all():
                raise ValueError('Expected nonnegative integer unique-read counts')
            joined = chunk[chunk.junction_id.isin(junctions)].merge(meta, left_on='sample_id', right_on='CB', validate='many_to_one')
            pieces.append(joined.groupby(keys).reads.sum())
        if not pieces:
            raise ValueError('Empty read table')
        counts = pd.concat(pieces).groupby(level=keys).sum()
        unit_keys = ['final_region', 'Group_L2', 'donor_id']
        observations = []
        for key, donor_meta in meta.groupby(unit_keys):
            try:
                donor = counts.loc[key]
            except KeyError:
                donor = pd.Series(dtype=float)
            for event in events:
                inc = np.mean([donor.get(j, 0) for j in event['inclusion']])
                exc = np.mean([donor.get(j, 0) for j in event['exclusion']])
                coverage = inc + exc
                if coverage < a.min_reads:
                    continue
                observations.append(dict(region=key[0], group=key[1], donor_id=key[2], event_id=event['event_id'], event_type=event['event_type'], gene_name=event['gene_name'], inclusion_support=inc, exclusion_support=exc, PSI=inc / coverage, n_cells=len(donor_meta), read_origin='original_STAR_unique'))
        obs = pd.DataFrame(observations)
        if obs.empty:
            raise ValueError('No donor/event meets the prespecified read threshold')
        obs.to_csv(out / 'donor_event_PSI.csv', index=False)
        rng = np.random.default_rng(a.seed)
        rows = []
        for (region, event), frame in obs.groupby(['region', 'event_id']):
            x = frame.loc[frame.group == 'DBS_I', 'PSI'].to_numpy()
            y = frame.loc[frame.group == 'Sham_I', 'PSI'].to_numpy()
            pv = np.nan
            if min(len(x), len(y)) >= a.min_donors:
                test = permutation_test((x, y), lambda a, b: np.mean(a) - np.mean(b), vectorized=False, permutation_type='independent', alternative='two-sided', n_resamples=a.permutations, random_state=rng)
                pv = test.pvalue
            mx = x.mean() if len(x) else np.nan
            my = y.mean() if len(y) else np.nan
            rows.append(dict(event_id=event, gene_name=frame.gene_name.iloc[0], event_type=frame.event_type.iloc[0], region=region, PSI_Sham=my, PSI_DBS=mx, delta_PSI=mx - my, P=pv, n_donors_Sham=len(y), n_donors_DBS=len(x), test='donor_permutation', read_origin='original_STAR_unique'))
        result = pd.DataFrame(rows)
        result['FDR'] = bh(result.P)
        result.to_csv(out / 'DAS_all_regions_donor.csv', index=False)
        summary = result.groupby('region').agg(n_events=('event_id', 'size'), n_tested=('P', 'count'))
        selected = result[(result.FDR < 0.05) & (result.delta_PSI.abs() >= 0.1)]
        summary['n_DAS'] = selected.groupby('region').size().reindex(summary.index, fill_value=0)
        summary.to_csv(out / 'region_DAS_counts_donor.csv')
    main()
    return locals()

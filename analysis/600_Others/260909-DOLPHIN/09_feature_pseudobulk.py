def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260909-DOLPHIN')
        import anndata as ad
        import numpy as np
        import pandas as pd
        from scipy import sparse
        from scipy.io import mmwrite
        x = inputs['x_4'][str(input_path(data, a.counts))]
        if not x.obs_names.is_unique or not x.var_names.is_unique or x.obs_names.str.contains('.aggr', regex=False).any():
            raise ValueError('Require unique original cell and feature IDs')
        if x.uns.get('read_origin') != 'original_unaggregated':
            raise ValueError('Require verified read_origin=original_unaggregated in H5AD.uns')
        counts = sparse.csr_matrix(x.X)
        if (counts.data < 0).any() or not np.equal(counts.data, np.round(counts.data)).all():
            raise ValueError('Expected original nonnegative integer counts')
        m = inputs['m_5'][str(input_path(data, a.metadata))].set_index('CB').loc[x.obs_names]
        keys = ['final_region', 'Group_L2', 'donor_id']
        if m[keys].isna().any().any() or m.donor_id.isin(['unassigned', 'doublet', '']).any():
            raise ValueError('Resolved biological donors are required')
        m['unit_id'] = m[keys].astype(str).agg('|'.join, axis=1)
        units = sorted(m.unit_id.unique())
        indicator = sparse.csr_matrix((np.ones(len(m)), (pd.Categorical(m.unit_id, categories=units).codes, np.arange(len(m)))), shape=(len(units), len(m)))
        aggregate = indicator @ counts
        target = out / (a.feature_kind + '_pseudobulk')
        target.mkdir(exist_ok=True)
        mmwrite(target / 'counts.mtx', aggregate.T)
        unit_meta = m.drop_duplicates('unit_id').set_index('unit_id').loc[units, keys]
        unit_meta.to_csv(target / 'samples.csv')
        var = x.var.copy()
        var.insert(0, 'feature_id', x.var_names)
        var.to_csv(target / 'features.csv', index=False)
    main()
    return locals()

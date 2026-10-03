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
        import pandas as pd
        from sklearn.neighbors import NearestNeighbors
        x = inputs['x_4'][str(input_path(data, a.latent))]
        m = inputs['m_5'][str(input_path(data, a.metadata))].set_index('CB')
        if not m.index.is_unique or not x.obs_names.is_unique:
            raise ValueError('Duplicate cell IDs')
        m = m.loc[x.obs_names]
        keys = ['donor_id', 'Group_L2', 'final_region']
        if m[keys].isna().any().any() or m.donor_id.isin(['unassigned', 'doublet', '']).any():
            raise ValueError('Resolved donor, condition, and region are required')
        rows = []
        for key, frame in m.groupby(keys, observed=True):
            index = x.obs_names.get_indexer(frame.index)
            z = x.obsm['X_z'][index]
            nn = NearestNeighbors(n_neighbors=min(a.neighbors + 1, len(z))).fit(z)
            distance, nearest = nn.kneighbors(z)
            for i, cell in enumerate(frame.index):
                rank = 0
                for dist, j in zip(distance[i], nearest[i]):
                    if j == i:
                        continue
                    rank += 1
                    rows.append(dict(CB=cell, neighbor=frame.index[j], rank=rank, distance=dist, donor_id=key[0], Group_L2=key[1], final_region=key[2]))
        pd.DataFrame(rows).to_csv(out / 'neighbors_within_donor.csv', index=False)
    main()
    return locals()

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    from pathlib import Path
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('500_GRN')
        import pandas as pd
        import numpy as np
        from scipy.io import mmread
        import loompy
        src = input_path(data, a.counts_dir) if a.counts_dir else out / 'counts'
        x = inputs['x_4'][str(src / 'counts.mtx')].tocsc()
        m = inputs['m_5'][str(src / 'metadata.csv')]
        genes = inputs['genes_6'][str(src / 'genes.csv')].iloc[:, 0].astype(str)
        if x.shape != (len(genes), len(m)) or genes.duplicated().any() or m.Cell_ID.duplicated().any():
            raise ValueError('Counts, unique genes, and unique cell IDs must align')
        if (x.data < 0).any() or not np.equal(x.data, np.round(x.data)).all():
            raise ValueError('Expected raw counts')
        attrs = {k: m[k].fillna('').to_numpy() for k in m.columns}
        attrs['CellID'] = m.Cell_ID.to_numpy()
        attrs['nGene'] = np.asarray((x > 0).sum(axis=0)).ravel()
        attrs['nUMI'] = np.asarray(x.sum(axis=0)).ravel()
        target = out / 'expression.loom'
        if inputs['data_7']:
            raise FileExistsError(target)
        loompy.create(str(target), x, {'Gene': genes.to_numpy()}, attrs)
    main()
    return locals()

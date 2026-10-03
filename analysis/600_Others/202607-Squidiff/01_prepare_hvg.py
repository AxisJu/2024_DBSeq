def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/202607-Squidiff')
        import anndata as ad
        import numpy as np
        import scanpy as sc
        from scipy import sparse
        x = inputs['x_4'][str(input_path(data, a.training))]
        x = x[x.obs.Group_L2.isin(['DBS_I', 'Sham_I'])].copy()
        x.obs['Group'] = x.obs.Group_L2.astype(str)
        if not x.var_names.is_unique or not x.obs_names.is_unique:
            raise ValueError('Gene and cell IDs must be unique')
        if a.frozen_features:
            ref = inputs['ref_5'][str(input_path(data, a.frozen_features))]
            genes = ref.var_names.tolist()
            if not set(genes) <= set(x.var_names):
                raise ValueError('Frozen checkpoint features are absent from training data')
        else:
            layer = None if a.counts_layer == 'X' else a.counts_layer
            counts = x.layers[layer] if layer else x.X
            values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
            if not np.isfinite(values).all() or (values < 0).any() or (not np.allclose(values, np.round(values))):
                raise ValueError('seurat_v3 requires raw counts; supply --counts-layer or --frozen-features')
            sc.pp.highly_variable_genes(x, n_top_genes=a.n_genes, flavor='seurat_v3', layer=layer)
            genes = x.var_names[x.var.highly_variable].tolist()
        x[:, genes].copy().write_h5ad(out / 'training_hvg.h5ad')
        import pandas as pd
        pd.DataFrame({'gene': genes}).to_csv(out / 'training_genes.csv', index=False)
    main()
    return locals()

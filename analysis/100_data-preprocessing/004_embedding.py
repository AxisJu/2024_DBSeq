def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']

    def main():
        p = parser(__doc__)
        args = inputs['args_3']
        paths = Paths(args)
        import os
        import scvi
        import torch
        import numpy as np
        import scanpy as sc
        import pandas as pd
        import anndata as ad
        scvi.settings.seed = args.seed
        torch.set_float32_matmul_precision('high')
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi.h5ad'))]
        adata = adata[:, adata.var.highly_variable]
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)
        adata.layers['scale'] = adata.X.copy()
        sc.tl.pca(adata, layer='scVI', svd_solver='arpack')
        sc.pp.neighbors(adata, knn=True, n_neighbors=30, n_pcs=30, use_rep='X_scVI')
        sc.tl.umap(adata, min_dist=1, spread=1, n_components=2)
        adata.write(paths.output('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap.h5ad'), compression='gzip')
        matrix = adata.layers['scale']
        np.savez(paths.output('data/snRNAseq_mouse/processed/matrix/matrix_scale_hvg.npz'), matrix)
        mat = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)
        mat.to_csv(paths.output('data/snRNAseq_mouse/processed/metadata/reduction_PCA.csv'), index=True)
        mat = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names)
        mat.to_csv(paths.output('data/snRNAseq_mouse/processed/metadata/reduction_UMAP_scVI.csv'), index=True)
    main()
    return locals()

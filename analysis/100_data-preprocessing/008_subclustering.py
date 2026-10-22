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
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L1.h5ad'))]
        adata_ExN = adata[adata.obs['celltype_level1'].isin(['ExN'])].copy()
        sc.pp.neighbors(adata_ExN, knn=True, n_neighbors=30, n_pcs=30, use_rep='X_scVI')
        sc.tl.umap(adata_ExN, min_dist=1, spread=1, n_components=2)
        for resolution in [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]:
            sc.tl.leiden(adata_ExN, resolution=resolution, key_added=f'leiden_ExN_res_{str(resolution)}')
        sc.tl.rank_genes_groups(adata_ExN, groupby='leiden_ExN_res_0.3', method='wilcoxon')
        marker_genes_df = sc.get.rank_genes_groups_df(adata_ExN, group=None)
        marker_genes_df.to_csv(paths.output('results/Table/markers_ExN_leidenres0.3_250624.csv'), index=False)
        adata_ExN.write(paths.output('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L1-ExN.h5ad'), compression='gzip')
        adata_InN = adata[adata.obs['celltype_level1'].isin(['InN'])].copy()
        sc.pp.neighbors(adata_InN, knn=True, n_neighbors=30, n_pcs=30, use_rep='X_scVI')
        sc.tl.umap(adata_InN, min_dist=1, spread=1, n_components=2)
        for resolution in [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]:
            sc.tl.leiden(adata_InN, resolution=resolution, key_added=f'leiden_InN_res_{str(resolution)}')
        sc.tl.rank_genes_groups(adata_InN, groupby='leiden_InN_res_0.4', method='wilcoxon')
        marker_genes_df = sc.get.rank_genes_groups_df(adata_InN, group=None)
        marker_genes_df.to_csv(paths.output('results/Table/markers_InN_leidenres0.4_250624.csv'), index=False)
        adata_InN.write(paths.output('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L1-InN.h5ad'), compression='gzip')
    main()
    return locals()

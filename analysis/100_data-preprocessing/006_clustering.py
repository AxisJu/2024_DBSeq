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
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap.h5ad'))]
        for resolution in [0.1, 0.3, 0.5, 0.7, 1.0, 2.0]:
            sc.tl.leiden(adata, resolution=resolution, key_added=f'leiden_res_{str(resolution)}')
        sc.tl.rank_genes_groups(adata, groupby='leiden_res_0.3', layer='scVI', method='wilcoxon')
        marker_genes_df = sc.get.rank_genes_groups_df(adata, group=None)
        marker_genes_df.to_csv(paths.output('results/250226 Clustering/markers_leidenres0.3_250217.csv'), index=False)
        adata.write(paths.output('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden.h5ad'), compression='gzip')
    main()
    return locals()

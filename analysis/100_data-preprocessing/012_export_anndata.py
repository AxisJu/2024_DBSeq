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
        import scipy.sparse as sp
        scvi.settings.seed = args.seed
        torch.set_float32_matmul_precision('high')
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_scanvi.h5ad'))]
        adata = adata.raw.to_adata()
        df = inputs['df_5'][str(paths.input('data/snRNAseq_mouse/processed/metadata/metadata_demultiplexing.csv'))]
        df = df.loc[adata.obs.index]
        for col in df.columns:
            adata.obs[col] = df[col].values
        counts = inputs['counts_6'][str(paths.input('data/snRNAseq_mouse/processed/matrix/matrix_counts.npz'))].transpose()
        adata.layers['counts'] = counts
        logdata = inputs['logdata_7'][str(paths.input('data/snRNAseq_mouse/processed/matrix/matrix_logCPM.npz'))].transpose()
        adata.layers['logCPM'] = logdata
        adata.write(paths.output('data/snRNAseq_mouse/processed/matrix/running_all_250704.h5ad'), compression='gzip')
    main()
    return locals()

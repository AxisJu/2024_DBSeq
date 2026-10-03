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
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden.h5ad'))]
        adata.obs['celltype_level1'] = adata.obs['leiden_res_0.3'].map({'1': 'ExN', '6': 'ExN', '10': 'ExN', '12': 'ExN', '13': 'ExN', '14': 'ExN', '16': 'ExN', '17': 'ExN', '18': 'ExN', '19': 'ExN', '20': 'ExN', '29': 'ExN', '25': 'ExN', '2': 'InN', '3': 'InN', '8': 'InN', '9': 'InN', '11': 'InN', '21': 'InN', '22': 'InN', '23': 'InN', '24': 'InN', '28': 'InN', '7': 'OPC', '0': 'Oligo', '4': 'Micro', '5': 'Astro', '15': 'VCs', '26': 'EPCs', '27': 'CHPCs'})
        adata.obs['celltype_level1'] = adata.obs['celltype_level1'].astype('category')
        adata.write(paths.output('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L1.h5ad'), compression='gzip')
        mat = pd.DataFrame(adata.obs['celltype_level1'], index=adata.obs_names)
        mat.to_csv(paths.output('data/snRNAseq_mouse/processed/metadata/metadata_annotation_celltype_level1.csv'), index=True)
    main()
    return locals()

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
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed.h5ad'))]
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True, layer='counts', flavor='seurat_v3', batch_key='Sample')
        scvi.model.SCVI.setup_anndata(adata, layer='counts', categorical_covariate_keys=['Sample_Time', 'Sequencing_Time'], continuous_covariate_keys=['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribo'])
        model = scvi.model.SCVI(adata, dropout_rate=0.0005, n_hidden=128, n_latent=30, n_layers=2, gene_likelihood='zinb')
        model.train(accelerator=args.accelerator, devices=1, max_epochs=1000, early_stopping=True)
        model.save(paths.output('data/snRNAseq_mouse/processed/intermediate/scVI/model_scVI_3000HVG2Layer30Latent126Hidden_250227'), overwrite=True)
        adata.obsm['X_scVI'] = model.get_latent_representation()
        adata.layers['scVI'] = model.get_normalized_expression(library_size=100000.0)
        adata.write(paths.output('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi.h5ad'), compression='gzip')
        matrix = adata.layers['scVI']
        np.savez(paths.output('data/snRNAseq_mouse/processed/matrix/matrix_scVI_hvg.npz'), matrix)
        mat = pd.DataFrame(adata.obsm['X_scVI'], index=adata.obs_names)
        mat.to_csv(paths.output('data/snRNAseq_mouse/processed/metadata/reduction_scVI.csv'), index=True)
    main()
    return locals()

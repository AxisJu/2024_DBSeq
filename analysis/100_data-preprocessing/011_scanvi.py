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
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L23.h5ad'))]
        model = inputs['model_5'][str(paths.input('data/snRNAseq_mouse/processed/intermediate/scVI/model_scVI_3000HVG2Layer30Latent126Hidden_250227'))]
        scanvi_model = scvi.model.SCANVI.from_scvi_model(model, adata=adata, labels_key='celltype_level2', unlabeled_category='Unknown')
        scanvi_model.train(accelerator=args.accelerator, devices=1, max_epochs=30, n_samples_per_label=1000)
        scanvi_model.save(paths.output('data/snRNAseq_mouse/processed/intermediate/scVI/model_scANVI_3000HVG_celltypeL2_250702'), overwrite=True)
        adata.obsm['X_scANVI'] = scanvi_model.get_latent_representation(adata)
        adata.layers['scANVI'] = scanvi_model.get_normalized_expression(library_size=100000.0)
        sc.pp.neighbors(adata, knn=True, n_neighbors=30, n_pcs=30, use_rep='X_scANVI')
        sc.tl.paga(adata, groups='celltype_level1')
        from scanpy.plotting._tools.paga import _compute_pos
        adjacency = adata.uns['paga']['connectivities'].copy()
        adjacency.data[adjacency.data < 0.01] = 0
        adjacency.eliminate_zeros()
        adata.uns['paga']['pos'] = _compute_pos(adjacency, layout='fr', random_state=0)
        sc.tl.umap(adata, min_dist=0.95, spread=0.95, n_components=2, init_pos='paga')
        adata.write(paths.output('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_scanvi.h5ad'), compression='gzip')
        adata.obs.to_csv(paths.output('data/snRNAseq_mouse/processed/metadata/metadata_obs_109.csv'), index=True)
        mat = pd.DataFrame(adata.obsm['X_scANVI'], index=adata.obs_names)
        mat.to_csv(paths.output('data/snRNAseq_mouse/processed/metadata/reduction_scANVI.csv'), index=True)
        mat = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names)
        mat.to_csv(paths.output('data/snRNAseq_mouse/processed/metadata/reduction_UMAP_scANVI.csv'), index=True)
    main()
    return locals()

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
        import anndata
        import scanpy as sc
        import pandas as pd
        sc.settings.verbosity = 1
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed.h5ad'))]
        sc.tl.pca(adata)
        sc.external.pp.bbknn(adata, batch_key='Sample', use_rep='X_pca')
        sc.tl.umap(adata, min_dist=1, spread=1, n_components=2)
        adata.uns['neighbors_BBKNN'] = adata.uns.pop('neighbors')
        adata.uns['umap_BBKNN'] = adata.uns.pop('umap')
        adata.obsm['X_umap_BBKNN'] = adata.obsm.pop('X_umap')
        adata.obsp['distances_BBKNN'] = adata.obsp.pop('distances')
        adata.obsp['connectivities_BBKNN'] = adata.obsp.pop('connectivities')
        adata.write(paths.output('data/snRNAseq_mouse/processed/scanpy/merged_preprocessed_bbknn.h5ad'), compression='gzip')
    main()
    return locals()

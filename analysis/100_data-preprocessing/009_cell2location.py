def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']
    MAIN_REGION = {'OPC': '/', 'Oligo': '/', 'Micro': '/', 'Astro': '/', 'VCs': '/', 'EPCs': '/', 'CHPCs': '/', 'ExN_0': '/', 'ExN_1': 'ATN', 'ExN_2': '/', 'ExN_3': 'RSP', 'ExN_4': 'RSP', 'ExN_5': 'HPC', 'ExN_6': '/', 'ExN_7': '/', 'ExN_8': '/', 'ExN_9': 'ATN', 'ExN_10': 'HPC', 'ExN_11': '/', 'ExN_12': 'HPC', 'ExN_13': '/', 'ExN_14': 'ATN', 'ExN_15': 'RSP', 'ExN_16': '/', 'ExN_17': '/', 'ExN_18': '/', 'InN_0': '/', 'InN_1': '/', 'InN_2': '/', 'InN_3': '/', 'InN_4': '/', 'InN_5': '/', 'InN_6': '/', 'InN_7': '/', 'InN_8': '/', 'InN_9': 'RSP', 'InN_10': '/', 'InN_11': '/', 'InN_12': 'ATN', 'InN_13': '/', 'InN_14': '/', 'InN_15': '/', 'InN_16': '/', 'InN_17': '/'}

    def spatial_prep(args, paths):
        import numpy as np
        import pandas as pd
        import anndata
        from scipy import sparse
        from sklearn.cluster import MiniBatchKMeans
        obj = inputs['obj_3'][str(paths.input(args.input))]
        annotation = inputs['annotation_4'][str(paths.input(args.annotation))]
        annotation['cell_uid'] = annotation.section_id.astype(str).str.replace('^T', '', regex=True) + '_' + annotation.cell_id.astype(str)
        annotation = annotation.set_index('cell_uid')
        if not annotation.index.is_unique:
            raise ValueError('Spatial cell_uid values are duplicated')
        ids = args.section.removeprefix('T') + '_' + obj.obs_names.astype(str)
        for key in ('cell_class', 'cell_subclass', 'cell_cluster'):
            obj.obs[key] = annotation[key].reindex(ids).to_numpy()
        obj = obj[obj.obs.cell_class.str.contains('Neurons', na=False)].copy()
        if obj.n_obs < args.cells_per_bin:
            raise ValueError('Insufficient neurons for spatial binning')
        model = MiniBatchKMeans(n_clusters=obj.n_obs // args.cells_per_bin, batch_size=10000, random_state=0, n_init=3)
        obj.obs['super_id'] = model.fit_predict(obj.obsm['spatial'])
        labels = np.sort(obj.obs.super_id.unique())
        incidence = sparse.csr_matrix((np.ones(obj.n_obs), (obj.obs.super_id.to_numpy(), np.arange(obj.n_obs))), shape=(len(labels), obj.n_obs))
        matrix = incidence @ sparse.csr_matrix(obj.X)

        def mode(values):
            count = values.value_counts()
            return count.index[0] if len(count) else 'unassigned'
        aggregations = {key: 'mean' for key in ('rx', 'ry')}
        aggregations.update({key: mode for key in ('gene_area', 'area_name', 'cell_class', 'cell_subclass', 'cell_cluster')})
        obs = obj.obs.groupby('super_id', sort=True).agg(aggregations).reindex(labels)
        obs.index = [f'super_{x}' for x in labels]
        result = anndata.AnnData(matrix, obs=obs, var=obj.var.copy())
        result.obsm['spatial'] = obs[['rx', 'ry']].to_numpy()
        mito = result.var_names.str.startswith('mt-')
        result.obsm['MT'] = result[:, mito].X.toarray()
        result = result[:, ~mito].copy()
        result.write_h5ad(paths.output(f'data/snRNAseq_mouse/processed/intermediate/cell2loaction/{args.section}_10Cellbin.h5ad'))

    def reference(args, paths):
        import numpy as np
        import anndata
        import cell2location
        from cell2location.utils.filtering import filter_genes
        labels = inputs['labels_5'][str(paths.input(args.input))]
        counts = inputs['counts_6'][str(paths.input(args.counts_input))]
        if not labels.obs_names.isin(counts.obs_names).all():
            raise ValueError('Reference labels cannot be matched to raw count cells')
        obj = counts[labels.obs_names].copy()
        obj.obs = labels.obs.copy()
        if args.lineage:
            obj = obj[obj.obs.celltype_level1 == args.lineage].copy()
        if args.group:
            obj = obj[obj.obs.Group_L2.isin(args.group)].copy()
        if args.label == 'mainregion':
            obj.obs['mainregion'] = obj.obs.leiden_celltype_level2.astype(str).map(MAIN_REGION).astype('category')
        if obj.obs[args.label].isna().any():
            raise ValueError('Reference labels contain missing values')
        if 'counts' in obj.layers:
            obj.X = obj.layers['counts'].copy()
        values = obj.X.data if hasattr(obj.X, 'indptr') else obj.X
        if not np.allclose(values, np.round(values)):
            raise ValueError('Reference requires raw counts, not normalized expression')
        selected = filter_genes(obj, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
        obj = obj[:, selected].copy()
        cell2location.models.RegressionModel.setup_anndata(obj, labels_key=args.label, categorical_covariate_keys=['Sample_Time', 'Sequencing_Time'], continuous_covariate_keys=['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribo'])
        model = cell2location.models.RegressionModel(obj)
        model.train(max_epochs=args.epochs or 300, accelerator=args.accelerator, device=args.device)
        obj = model.export_posterior(obj, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': args.accelerator != 'cpu'})
        base = 'data/snRNAseq_mouse/processed/intermediate/cell2loaction/reference_signatures'
        model.save(str(paths.output(f'{base}/model_{args.name}')), overwrite=args.force)
        obj.write_h5ad(paths.output(f'{base}/sc_{args.name}.h5ad'))

    def mapping(args, paths):
        import numpy as np
        import anndata
        import cell2location
        ref = inputs['ref_7'][str(paths.input(args.reference))]
        obj = inputs['obj_8'][str(paths.input(args.input))]
        names = ref.uns['mod']['factor_names']
        columns = [f'means_per_cluster_mu_fg_{name}' for name in names]
        signatures = ref.varm['means_per_cluster_mu_fg'][columns].copy() if 'means_per_cluster_mu_fg' in ref.varm else ref.var[columns].copy()
        signatures.columns = names
        common = np.intersect1d(obj.var_names, signatures.index)
        if not len(common):
            raise ValueError('Spatial and reference objects have no shared genes')
        obj = obj[:, common].copy()
        signatures = signatures.loc[common]
        cell2location.models.Cell2location.setup_anndata(obj)
        model = cell2location.models.Cell2location(obj, cell_state_df=signatures, N_cells_per_location=args.cells_per_bin, detection_alpha=20)
        model.train(max_epochs=args.epochs or 15000, batch_size=None, train_size=1, accelerator=args.accelerator, device=args.device)
        obj = model.export_posterior(obj, sample_kwargs={'num_samples': 1000, 'batch_size': obj.n_obs, 'use_gpu': args.accelerator != 'cpu'})
        base = 'data/snRNAseq_mouse/processed/intermediate/cell2loaction/cell2location_map'
        model.save(str(paths.output(f'{base}/model_{args.name}_{args.section}')), overwrite=args.force)
        obj.write_h5ad(paths.output(f'{base}/sp_{args.name}_{args.section}.h5ad'))

    def main():
        p = parser(__doc__)
        args = inputs['args_9']
        paths = Paths(args)
        if args.stage == 'spatial-prep' and (not args.annotation):
            p.error('spatial-prep requires --annotation')
        if args.stage == 'map' and (not args.reference):
            p.error('map requires --reference')
        import numpy as np
        np.random.seed(args.seed)
        if args.stage == 'spatial-prep':
            spatial_prep(args, paths)
        else:
            import scvi
            scvi.settings.seed = args.seed
            {'reference': reference, 'map': mapping}[args.stage](args, paths)
    main()
    return locals()

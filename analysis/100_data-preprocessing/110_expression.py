def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']

    def library_scaled_log(counts, samples):
        import numpy as np
        from scipy import sparse
        x = sparse.csc_matrix(counts, dtype=float)
        totals = np.asarray(x.sum(axis=0)).ravel()
        if x.shape[1] != len(samples) or np.any(totals <= 0):
            raise ValueError('Each metadata row must match a nonempty count column')
        samples = np.asarray(samples)
        targets = np.zeros(len(samples), dtype=float)
        for sample in np.unique(samples):
            mask = samples == sample
            targets[mask] = np.median(totals[mask])
        x = x @ sparse.diags(targets / totals)
        x.data = np.log2(x.data + 1)
        return x.tocsc()

    def main():
        p = parser(__doc__)
        args = inputs['args_3']
        paths = Paths(args)
        import numpy as np
        import pandas as pd
        from scipy import sparse
        base = 'data/snRNAseq_mouse/processed'
        if args.method == 'logcpm':
            counts = inputs['counts_4'][str(paths.input(base + '/matrix/matrix_counts.npz'))]
            meta = inputs['meta_5'][str(paths.input(base + '/metadata/metadata_obs_101.csv'))]
            result = library_scaled_log(counts, meta.Sample)
            sparse.save_npz(paths.output(base + '/matrix/matrix_logCPM.npz'), result)
        else:
            import anndata
            obj = inputs['obj_6'][str(paths.input(base + '/matrix/merged_preprocessed.h5ad'))]
            if args.method == 'magic':
                import magic
                import scprep
                values = scprep.transform.sqrt(scprep.normalize.library_size_normalize(obj.layers['counts']))
                model = magic.MAGIC(knn=5, knn_max=15, decay=1, t=3, random_state=args.seed)
                np.savez(paths.output(base + '/matrix/matrix_MAGIC.npz'), model.fit_transform(values, genes='all_genes'))
            else:
                import os
                import scNET
                scNET.main.MAX_CELLS_BATCH_SIZE = 3000
                scNET.main.NUM_LAYERS = 5
                scNET.main.DE_GENES_NUM = 3500
                work = paths.output(base + '/intermediate/scNET')
                work.mkdir(parents=True, exist_ok=True)
                os.chdir(work)
                if args.train_scnet:
                    scNET.run_scNET(obj, pre_processing_flag=False, human_flag=False, number_of_batches=20, split_cells=False, max_epoch=300, model_name=args.model_name, save_model_flag=True)
                _, _, nodes, outputs = inputs['nodes_outputs_7'][str(args.model_name)]
                reconstructed = scNET.create_reconstructed_obj(nodes, outputs, obj)
                reconstructed.write_h5ad(paths.output(base + '/intermediate/scNET/matrix_scNET.h5ad'), compression='gzip')
                sparse.save_npz(paths.output(base + '/matrix/matrix_scNET.npz'), sparse.csc_matrix(reconstructed.X.T))
                reconstructed.obs.to_csv(paths.output(base + '/metadata/metadata_obs_scNET.csv'))
                reconstructed.var.to_csv(paths.output(base + '/metadata/metadata_var_scNET.csv'))
    main()
    return locals()

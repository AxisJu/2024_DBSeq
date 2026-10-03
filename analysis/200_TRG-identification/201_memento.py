def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']

    def select_contrast(adata, column, celltype, treatment, reference):
        mask = (adata.obs[column].astype(str) == str(celltype)) & adata.obs.Group_L2.isin([treatment, reference])
        selected = adata[mask].copy()
        selected.obs['Group'] = (selected.obs.Group_L2 == treatment).astype(int)
        if set(selected.obs.Group_L2.unique()) != {treatment, reference}:
            raise ValueError('A contrast must contain exactly the treatment and reference groups')
        return selected

    def main():
        p = parser(__doc__)
        args = inputs['args_3']
        paths = Paths(args)
        import anndata as ad
        import numpy as np
        import pandas as pd
        import memento
        from scipy import sparse
        from statsmodels.stats.multitest import multipletests
        np.random.seed(args.seed)
        if args.scope == 'atn':
            base = 'data/snRNAseq_mouse/processed'
            exp = inputs['exp_4'][str(paths.input(base + '/matrix/running_ATN_250704_mat.csv'))].T
            meta = inputs['meta_5'][str(paths.input(base + '/metadata/running_ATN_250704_metadata.csv'))]
            if not exp.index.is_unique or not meta.index.is_unique or (not exp.index.isin(meta.index).all()):
                raise ValueError('ATN expression and metadata cell IDs do not align')
            obj = ad.AnnData(sparse.csr_matrix(exp.values), obs=meta.loc[exp.index], var=pd.DataFrame(index=exp.columns))
            obj.layers['counts'] = obj.X.copy()
            columns = {'ANT': 'celltype'}
        else:
            obj = inputs['obj_6'][str(paths.input(args.input))]
            columns = {f'Level{level}': f'celltype_level{level}' for level in args.levels}
        comparisons = {'DBSvSham': ('DBS_I', 'Sham_I'), 'ShamvSaline': ('Sham_I', 'Saline_I'), 'DBSIvC': ('DBS_I', 'DBS_C'), 'ShamIvPTZ': ('Sham_I', 'PTZ')}
        if args.scope == 'atn':
            comparisons = {name: comparisons[name] for name in ('DBSvSham', 'ShamvSaline')}
        minimum = args.min_cells if args.min_cells is not None else 1 if args.scope == 'atn' else 50
        audit = []
        for level, column in columns.items():
            for celltype in sorted(obj.obs[column].dropna().astype(str).unique()):
                for name, (treatment, reference) in comparisons.items():
                    counts = obj.obs.loc[obj.obs[column].astype(str) == celltype, 'Group_L2'].value_counts()
                    if min(counts.get(treatment, 0), counts.get(reference, 0)) <= minimum:
                        continue
                    selected = select_contrast(obj, column, celltype, treatment, reference)
                    if args.scope != 'atn':
                        keep = np.ones(selected.n_vars, dtype=bool)
                        for group in (reference, treatment):
                            values = selected[selected.obs.Group_L2 == group].X
                            keep &= np.asarray((values > 0).mean(axis=0)).ravel() > args.fraction
                        selected = selected[:, keep].copy()
                    selected.X = selected.layers['counts'].copy()
                    result = memento.binary_test_1d(adata=selected, capture_rate=0.15, treatment_col='Group', num_cpus=args.threads, num_boot=args.bootstrap, verbose=False)
                    result = result.sort_values('de_coef', ascending=False)
                    finite = np.isfinite(result.de_pval)
                    result['FDR'] = np.nan
                    result.loc[finite, 'FDR'] = multipletests(result.loc[finite, 'de_pval'], method='fdr_bh')[1]
                    root = 'results/Table/MEMENTO/ANT' if args.scope == 'atn' else f'results/Table/MEMENTO_v2/{level}'
                    result.to_csv(paths.output(f'{root}/MEMENTO_1d_{name}_{celltype}.csv'), index=False)
                    audit.append(dict(level=level, celltype=celltype, contrast=name, reference=reference, treatment=treatment, reference_cells=counts[reference], treatment_cells=counts[treatment], inference_unit='cell; donor not included in binary_test_1d'))
        pd.DataFrame(audit).to_csv(paths.output(f'results/Table/MEMENTO_{args.scope}_units.csv'), index=False)
    main()
    return locals()

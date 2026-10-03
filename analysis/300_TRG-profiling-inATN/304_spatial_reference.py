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
        import pandas as pd
        import subprocess
        from scipy import sparse
        base = 'data/snRNAseq_mouse/processed'
        expression = inputs['expression_4'][str(paths.input(base + '/matrix/running_ATN_250704_mat.csv'))].T
        metadata = inputs['metadata_5'][str(paths.input(base + '/metadata/running_ATN_250704_metadata.csv'))]
        if not expression.index.is_unique or not expression.columns.is_unique or (not expression.index.isin(metadata.index).all()):
            raise ValueError('ATN expression identifiers do not align with metadata')
        obj = anndata.AnnData(sparse.csr_matrix(expression.to_numpy()), obs=metadata.loc[expression.index], var=pd.DataFrame(index=expression.columns))
        obj.layers['counts'] = obj.X.copy()
        filename = base + '/matrix/running_ATN_spatial_reference.h5ad'
        obj.write_h5ad(paths.output(filename))
        if not args.prepare_only:
            script = Path(__file__).resolve().parent.parent / '100_data-preprocessing' / '108_cell2location.py'
            command = [sys.executable, str(script), 'reference', '--data-root', args.data_root, '--output-root', args.output_root, '--input', filename, '--counts-input', filename, '--label', 'celltype', '--name', 'ant-celltype', '--epochs', str(args.epochs), '--accelerator', args.accelerator, '--seed', str(args.seed), '--threads', str(args.threads)]
            if args.force:
                command.append('--force')
            inputs['data_6']
    main()
    return locals()

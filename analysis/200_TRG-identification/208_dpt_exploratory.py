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
        import scanpy as sc
        import numpy as np
        obj = inputs['obj_4'][str(paths.input(args.input))]
        obj = obj[obj.obs.celltype_level3.isin(args.celltypes) & obj.obs.Group_L2.isin(['DBS_I', 'Sham_I', 'Saline_I'])].copy()
        sc.tl.diffmap(obj)
        sc.pp.neighbors(obj, n_neighbors=args.neighbors, use_rep='X_diffmap', random_state=args.seed)
        values = obj.obsm['X_diffmap'][:, args.root_component]
        obj.uns['iroot'] = int(np.argmin(values) if args.root_extreme == 'min' else np.argmax(values))
        sc.tl.dpt(obj)
        obj.obs[['dpt_pseudotime']].to_csv(paths.output('results/Table/exploratory_dpt.csv'))
        obj.write_h5ad(paths.output('results/models/exploratory_dpt.h5ad'))
    main()
    return locals()

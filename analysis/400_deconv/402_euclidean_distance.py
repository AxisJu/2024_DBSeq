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
        import pandas as pd
        from scipy.spatial.distance import cdist
        raw = inputs['raw_4'][str(paths.input(args.coordinates))]
        coords = raw[['sc213_acr', 'x', 'y', 'z']].dropna().groupby('sc213_acr', sort=True)[['x', 'y', 'z']].mean()
        selected = coords.loc[args.targets]
        all_distances = pd.DataFrame(cdist(coords, coords), index=coords.index, columns=coords.index)
        target_distances = pd.DataFrame(cdist(coords, selected), index=coords.index, columns=[f'ED_to_{s}' for s in args.targets])
        base = 'data/snRNAseq_mouse/processed/intermediate/BCT'
        all_distances.to_csv(paths.output(base + '/ED_all.csv'), float_format='%.6f')
        target_distances.to_csv(paths.output(base + '/ED_ant.csv'), float_format='%.6f')
    main()
    return locals()

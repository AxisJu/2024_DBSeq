def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']

    def mean_connectivity(series):
        import numpy as np
        correlations = [np.corrcoef(np.asarray(values, dtype=float)) for values in series]
        if not correlations or len({matrix.shape for matrix in correlations}) != 1:
            raise ValueError('Each animal must have the same number of region rows')
        return np.mean(np.stack(correlations, axis=2), axis=2)

    def main():
        p = parser(__doc__)
        args = inputs['args_3']
        paths = Paths(args)
        import numpy as np
        from scipy.io import loadmat
        obj = inputs['obj_4'][str(paths.input(args.timeseries))][args.variable]
        matrices = [value for value in obj.ravel()]
        result = mean_connectivity(matrices)
        np.savetxt(paths.output('data/derivatives/dbseq/mouse_BOLD_fc.csv'), result, delimiter=',')
    main()
    return locals()

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import subprocess
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('500_GRN')
        exp = input_path(data, a.expression) if a.expression else out / 'expression.loom'
        tf, ranks, motifs = [input_path(data, s) for s in (a.tf_list, a.rankings, a.motifs)]
        for path in (exp, tf, ranks, motifs):
            if not inputs['data_4']:
                raise FileNotFoundError(path)
        common = ['--num_workers', str(a.workers)]
        commands = {'grn': [a.executable, 'grn', str(exp), str(tf), '--method', 'grnboost2', '--sparse', '--seed', str(a.seed), '-o', str(out / 'adj.tsv'), *common], 'ctx': [a.executable, 'ctx', str(out / 'adj.tsv'), str(ranks), '--annotations_fname', str(motifs), '--expression_mtx_fname', str(exp), '--mode', 'dask_multiprocessing', '--mask_dropouts', '-o', str(out / 'reg.csv'), *common], 'aucell': [a.executable, 'aucell', str(exp), str(out / 'reg.csv'), '--seed', str(a.seed), '-o', str(out / 'pyscenic_output.loom'), *common]}
        for stage in commands:
            if a.stage in ('all', stage):
                inputs['data_5']
    main()
    return locals()

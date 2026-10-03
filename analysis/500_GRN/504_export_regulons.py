def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import ast
    import csv
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def parse_targets(path):
        handle = inputs['handle_3'][str(path)]
        rows = csv.reader(handle)
        header = next(rows)
        while 'TargetGenes' not in header:
            header = next(rows)
        index = header.index('TargetGenes')
        for row in rows:
            if len(row) <= index or not row[index].startswith('['):
                continue
            for gene, weight in ast.literal_eval(row[index]):
                yield {'TF': row[0], 'motif': row[1], 'target_gene': str(gene), 'importance': float(weight)}

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_4']
        data, out = roots('500_GRN')
        import h5py
        import pandas as pd
        import numpy as np
        loom = input_path(data, a.loom) if a.loom else out / 'pyscenic_output.loom'
        reg = input_path(data, a.regulons) if a.regulons else out / 'reg.csv'
        f = inputs['f_5'][str(loom)]
        ca = f['col_attrs']
        decode = lambda x: [v.decode() if isinstance(v, bytes) else str(v) for v in x]
        ids = decode(ca['CellID'][:])
        raw = ca['RegulonsAUC'][:]
        if raw.dtype.names is None:
            raise ValueError('Expected named pySCENIC compound RegulonsAUC array')
        auc = pd.DataFrame({name: np.asarray(raw[name]).ravel() for name in raw.dtype.names}, index=ids)
        auc.index.name = 'Cell_ID'
        auc.to_csv(out / 'regulon_auc.csv')
        fields = ['donor_id', 'Group_L2', 'Sample', 'celltype_level1', 'final_region']
        metadata = {'Cell_ID': ids}
        for field in fields:
            if field in ca:
                metadata[field] = decode(ca[field][:])
        pd.DataFrame(metadata).to_csv(out / 'regulon_metadata.csv', index=False)
        pd.DataFrame(parse_targets(reg)).drop_duplicates().to_csv(out / 'regulon_targets.csv', index=False)
    main()
    return locals()

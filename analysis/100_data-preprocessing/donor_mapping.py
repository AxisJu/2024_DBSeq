def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    from pathlib import Path

    def read_cell_donors(data_root, output_root=None):
        import pandas as pd
        relative = Path('data/snRNAseq_mouse/processed/metadata')

        def locate(name):
            candidates = [Path(output_root) / relative / name] if output_root else []
            candidates.append(Path(data_root) / relative / name)
            return next((p for p in candidates if inputs['data_1']), candidates[-1])
        names = inputs['names_2'][str(locate('metadata_cellnamemapping.csv'))]
        obs = inputs['obs_3'][str(locate('metadata_obs_101.csv'))]
        donors = inputs['donors_4'][str(locate('metadata_demultiplexing.csv'))]
        for frame in (obs, donors):
            if not frame.index.is_unique:
                raise ValueError('Duplicate unique cell IDs in metadata')
        if names.unique_cellname.duplicated().any():
            raise ValueError('Duplicate unique_cellname in barcode mapping')
        merged = names.set_index('unique_cellname').join(obs[['Sample']], validate='one_to_one').join(donors[['donor_id']], validate='one_to_one')
        if merged[['Sample', 'donor_id']].isna().any().any():
            raise ValueError('Incomplete cell-to-donor join')
        if merged.duplicated(['Sample', 'raw_cellname']).any():
            raise ValueError('Sample and raw barcode must identify one cell')
        return merged.rename_axis('cell_id').reset_index()
    return locals()

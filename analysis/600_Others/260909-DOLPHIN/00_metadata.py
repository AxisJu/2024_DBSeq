def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260909-DOLPHIN')
        import pandas as pd
        base = input_path(data, a.metadata_dir)
        mapping = inputs['mapping_4'][str(base / 'metadata_cellnamemapping.csv')]
        demux = inputs['demux_5'][str(base / 'metadata_demultiplexing.csv')]
        columns = inputs['columns_6'][str(base / 'metadata_obs_101.csv')].columns
        obs = inputs['obs_7'][str(base / 'metadata_obs_101.csv')]
        joined = mapping.merge(demux[['donor_id']], left_on='unique_cellname', right_index=True, validate='one_to_one').merge(obs, left_on='unique_cellname', right_index=True, validate='one_to_one')
        joined['clean_barcode'] = joined.raw_cellname.str.extract('([ACGT]{16}-[0-9]+)$', expand=False)
        if joined.clean_barcode.isna().any() or joined.duplicated(['Sample', 'clean_barcode']).any():
            raise ValueError('Ambiguous sample/barcode mapping')
        joined['CB'] = joined.Sample.str.replace('_I_', '-', regex=False) + '_subset.TAG_CB_' + joined.clean_barcode
        joined = joined.rename(columns={'unique_cellname': 'Cell_ID'})
        meta = inputs['meta_8'][str(input_path(data, a.metadata))]
        meta = meta.drop(columns=['donor_id', 'Cell_ID', 'clean_barcode'], errors='ignore').merge(joined[['CB', 'Cell_ID', 'donor_id', 'clean_barcode']], on='CB', how='left', validate='one_to_one')
        if meta.donor_id.isna().any():
            raise ValueError('Some DOLPHIN cells do not match the donor mapping')
        excluded = meta.donor_id.isin(['unassigned', 'doublet', ''])
        meta.loc[excluded].to_csv(out / 'unresolved_donor_cells.csv', index=False)
        meta.loc[~excluded].to_csv(out / 'metadata_donor.csv', index=False)
    main()
    return locals()

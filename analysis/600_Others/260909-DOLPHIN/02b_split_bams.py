def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import subprocess
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260909-DOLPHIN')
        import pandas as pd
        samples = inputs['samples_4'][str(input_path(data, a.samples))]
        meta = inputs['meta_5'][str(input_path(data, a.metadata))]
        rows = []
        for sample in samples.itertuples():
            cells = meta[meta.Sample == sample.Sample]
            if cells.empty:
                continue
            path = out / 'split_bams' / sample.Sample
            path.mkdir(parents=True, exist_ok=True)
            barcode_file = path / 'barcodes.txt'
            cells.clean_barcode.to_csv(barcode_file, header=False, index=False)
            prefix = str(sample.Sample).replace('_I_', '-')
            subset = path / (prefix + '_subset.bam')
            inputs['data_6']
            inputs['data_7']
            for cell in cells.CB:
                bam = path / (cell + '.bam')
                if not inputs['data_8']:
                    raise FileNotFoundError(bam)
                rows.append({'CB': cell, 'bam': str(bam)})
        pd.DataFrame(rows).to_csv(out / 'cell_bam_manifest.csv', index=False)
    main()
    return locals()

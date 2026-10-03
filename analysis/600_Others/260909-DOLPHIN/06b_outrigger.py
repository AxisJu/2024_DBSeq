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
        rows = inputs['rows_4'][str(input_path(data, a.manifest))]
        if rows.sample_id.astype(str).str.contains('.aggr', regex=False).any():
            raise ValueError('Use original junctions for the inferential event catalog')
        paths = [str(input_path(data, value)) for value in rows.sj_path]
        target = out / 'outrigger_original'
        target.mkdir(exist_ok=True)
        inputs['data_5']
        inputs['data_6']
    main()
    return locals()

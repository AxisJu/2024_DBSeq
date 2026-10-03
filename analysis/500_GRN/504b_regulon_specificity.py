def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('500_GRN')
        import numpy as np
        import pandas as pd
        from scipy.spatial.distance import jensenshannon
        auc = inputs['auc_4'][str(input_path(data, a.auc) if a.auc else out / 'regulon_auc.csv')]
        meta = inputs['meta_5'][str(input_path(data, a.metadata) if a.metadata else out / 'regulon_metadata.csv')].loc[auc.index]
        if (auc.to_numpy() < 0).any() or not np.isfinite(auc.to_numpy()).all():
            raise ValueError('AUC must be finite and nonnegative')
        rows = []
        for group in sorted(meta[a.annotation].dropna().unique()):
            target = (meta[a.annotation] == group).to_numpy().astype(float)
            for regulon in auc.columns:
                value = auc[regulon].to_numpy()
                score = 1 - jensenshannon(value, target) if value.sum() > 0 else np.nan
                rows.append(dict(annotation=a.annotation, group=group, regulon=regulon, RSS=score))
        pd.DataFrame(rows).to_csv(out / ('regulon_specificity_' + a.annotation + '.csv'), index=False)
    main()
    return locals()

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import os
    from pathlib import Path
    import json
    import numpy as np
    import pandas as pd
    from scipy import stats
    import mpmath as mp
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def main():
        parser = argparse.ArgumentParser(description='Plot the frozen model calibration gene means.')
        args = inputs['args_1']
        path = args.input or Path(os.environ['DBSEQ_INPUT_PACKAGE_ROOT']) / 'sourcedata/figure 7/7d.csv'
        output = Path(os.environ['DBSEQ_OUTPUT_ROOT'])
        for kind in ['figures', 'sourcedata', 'statistics']:
            (output / kind / 'figure 7').mkdir(parents=True, exist_ok=True)
        data = inputs['data_2'][str(path)]
        assert not data.gene.duplicated().any()
        x = data.real_dbs_mean.to_numpy()
        y = data.virtual_dbs_scaled.to_numpy()
        assert np.isfinite(x).all() and np.isfinite(y).all()
        r = float(stats.pearsonr(x, y).statistic)
        n = len(x)
        mp.mp.dps = 60
        p = mp.betainc(mp.mpf(n - 2) / 2, mp.mpf('.5'), 0, 1 - mp.mpf(str(r)) ** 2, regularized=True)
        logp = float(mp.log10(p))
        evidence = dict(n_genes=n, Pearson_r=r, P=mp.nstr(p, 12), log10_P=logp, unit='gene means', interpretation='model calibration; not independent animal validation')
        (output / 'statistics/figure 7/7d_statistics.json').write_text(json.dumps(evidence, indent=2), encoding='utf-8')
        data.to_csv(output / 'sourcedata/figure 7/7d.csv', index=False)
        plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'pdf.fonttype': 42})
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
        ax.scatter(x, y, c='#ffd700', alpha=0.6, s=12, rasterized=True)
        limits = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(limits, limits, '--', color='grey', lw=0.8)
        ax.set(xlabel='Empirical ATN DBS mean expression', ylabel='Calibrated virtual DBS mean expression', xlim=limits, ylim=limits)
        ax.set_title(f'Pearson r={r:.3f}; log10(P)={logp:.2f}', fontsize=9)
        ax.set_aspect('equal')
        fig.tight_layout()
        for ext in ['png']:
            fig.savefig(output / f'figures/figure 7/7d.{ext}', dpi=400)
        print(json.dumps(evidence, indent=2))
    main()
    return locals()

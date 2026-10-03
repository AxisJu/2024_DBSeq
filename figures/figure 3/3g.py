def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import os
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy.spatial import cKDTree
    from scipy import stats
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    DATA = Path(os.environ['DBSEQ_DATA_ROOT'])
    OUTPUT = Path(os.environ['DBSEQ_OUTPUT_ROOT'])
    REGIONS = ['ANT', 'RSP', 'SUB']
    GROUPS = ['Sham', 'DBS']
    COLORS = ['#7e418f', '#ded33d']
    MIN_SNR = 3.0
    MIN_PLANES = 2
    ANCHOR_UM = 0.5
    NEIGHBOURHOOD_UM = 0.8

    def destination(kind, name):
        path = OUTPUT / kind / 'figure 3' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def observations():
        root = DATA / 'results/260813_spotiflow'
        inventory = inputs['inventory_1'][str(root / 'metadata/nd2_inventory.csv')]
        inventory = inventory[inventory.region.isin(REGIONS)].copy()
        assert not inventory.slide_id.duplicated().any()
        rows = []
        for item in inventory.itertuples(index=False):
            assert item.group in {'KA', 'DBS_KA'}
            path = root / 'spot_coordinates' / f'{item.slide_id}_spots.csv'
            spots = inputs['spots_2'][str(path)]
            qualified = spots[(spots.snr >= MIN_SNR) & (spots.n_slices >= MIN_PLANES)]
            channels = [qualified[qualified.channel_id.eq(c)] for c in [1, 2, 3]]
            psd, nr, nl = channels
            values = []
            if all((len(c) for c in channels)):
                xyz = ['x_um', 'y_um', 'z_um']
                distance, _ = cKDTree(psd[xyz]).query(nl[xyz])
                anchored = nl[distance <= ANCHOR_UM]
                neighbours = cKDTree(nr[xyz]).query_ball_point(anchored[xyz], NEIGHBOURHOOD_UM)
                for indexes, peak in zip(neighbours, anchored.intensity_peak):
                    if indexes:
                        values.append(float(peak * nr.iloc[indexes].intensity_peak.mean()))
            rows.append(dict(Slide_ID=item.slide_id, Mouse_ID=f'{item.group}_{item.mouse}', Brain_Region=item.region, Group='Sham' if item.group == 'KA' else 'DBS', Anchored_NLGN1_neighbourhoods=len(values), CoIntensity_1e6_au=np.mean(values) / 1000000.0 if values else np.nan, Status='observed' if values else 'no_qualified_neighbourhood', Source=str(path.relative_to(DATA))))
        return pd.DataFrame(rows)

    def main():
        raw = observations()
        raw.to_csv(destination('sourcedata', '3g_fov.csv'), index=False)
        animal = raw.groupby(['Brain_Region', 'Group', 'Mouse_ID'], as_index=False).agg(CoIntensity_1e6_au=('CoIntensity_1e6_au', 'mean'), FOV_with_metric=('CoIntensity_1e6_au', 'count'), FOV_total=('Slide_ID', 'size'))
        animal.to_csv(destination('sourcedata', '3g.csv'), index=False)
        results = []
        for region in REGIONS:
            sub = animal[animal.Brain_Region.eq(region)]
            sham, dbs = [sub.loc[sub.Group.eq(group), 'CoIntensity_1e6_au'].dropna().to_numpy() for group in GROUPS]
            p = stats.ttest_ind(dbs, sham, equal_var=False).pvalue if min(len(sham), len(dbs)) >= 2 else np.nan
            sd = np.sqrt(((len(dbs) - 1) * np.var(dbs, ddof=1) + (len(sham) - 1) * np.var(sham, ddof=1)) / (len(dbs) + len(sham) - 2))
            results.append(dict(Region=region, n_Sham=len(sham), n_DBS=len(dbs), P=p, CohenD=(dbs.mean() - sham.mean()) / sd if sd > 0 else np.nan))
        summary = pd.DataFrame(results)
        summary['FDR'] = stats.false_discovery_control(summary.P.to_numpy())
        summary.to_csv(destination('statistics', '3g_statistics.csv'), index=False)
        plt.rcParams.update({'font.family': 'Arial', 'pdf.fonttype': 42, 'font.size': 9})
        fig, ax = plt.subplots(figsize=(5, 4))
        for j, group in enumerate(GROUPS):
            for i, region in enumerate(REGIONS):
                values = animal.loc[animal.Group.eq(group) & animal.Brain_Region.eq(region), 'CoIntensity_1e6_au'].dropna().to_numpy()
                position = i + (j - 0.5) * 0.35
                ax.bar(position, values.mean(), width=0.32, color=COLORS[j], edgecolor='black', yerr=stats.sem(values), capsize=3, label=group if i == 0 else None)
                ax.scatter(position + np.linspace(-0.055, 0.055, len(values)), values, s=20, c='black' if j == 0 else 'white', edgecolors='black', zorder=3)
        ymax = float(animal.CoIntensity_1e6_au.max())
        for i, row in summary.iterrows():
            ax.text(i, ymax * 1.1, f'P={row.P:.3g}\nBH q={row.FDR:.3g}', ha='center', fontsize=8)
        ax.set(xticks=range(3), xticklabels=REGIONS, ylabel='Animal mean co-intensity (a.u. / 10$^6$)', ylim=(0, ymax * 1.3))
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False)
        fig.tight_layout()
        for ext in ['png']:
            fig.savefig(destination('figures', f'3g.{ext}'), dpi=400)
        print(summary.to_string(index=False))
    main()
    return locals()

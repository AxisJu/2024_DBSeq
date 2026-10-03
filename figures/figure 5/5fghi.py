def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import os
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    DATA = Path(os.environ['DBSEQ_DATA_ROOT'])
    OUTPUT = Path(os.environ['DBSEQ_OUTPUT_ROOT'])
    REGIONS = ['ATN', 'CP', 'RSPd', 'RT', 'RE', 'GPe', 'PF', 'DG', 'LSc', 'TRS']
    GROUPS = ['Sham', 'DBS']
    COLORS = ['#91459a', '#ffdc17']

    def destination(kind, name):
        path = OUTPUT / kind / 'figure 5' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def main():
        table = DATA / 'results/260814_Rsrc1-Sf3b1 IF/260929_Representative/V4/analysis/All_slices_NeuN_Target_Counts_Long_V4.csv'
        raw = inputs['raw_1'][str(table)]
        included = raw.region_present.astype(str).str.lower().eq('true') & raw.region.isin(REGIONS) & raw.group.isin(GROUPS)
        source = raw.loc[included, ['target', 'Sample_ID', 'mouse_id', 'group', 'region', 'n_NeuN_Target', 'ROI_area_mm2']].copy()
        assert not source.duplicated(['target', 'Sample_ID', 'region']).any()
        assert source.n_NeuN_Target.notna().all() and (source.n_NeuN_Target >= 0).all()
        animal = source.groupby(['target', 'mouse_id', 'group', 'region'], as_index=False).agg(Count=('n_NeuN_Target', 'mean'), n_sections=('Sample_ID', 'size'))
        sc2 = inputs['sc2_2'][str(DATA / 'results/260814_Rsrc1-Sf3b1 IF/260819_Sf3b1/results/tables/imaging_metrics_13_regions.csv')]
        plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'pdf.fonttype': 42})
        all_stats = []
        for target, panel_count, panel_corr in [('Sf3b1', 'f', 'g'), ('Rsrc1', 'h', 'i')]:
            frame = animal[animal.target.eq(target)]
            source[source.target.eq(target)].to_csv(destination('sourcedata', f'5{panel_count}_sections.csv'), index=False)
            frame.to_csv(destination('sourcedata', f'5{panel_count}.csv'), index=False)
            results = []
            for region in REGIONS:
                sub = frame[frame.region.eq(region)]
                sham, dbs = [sub.loc[sub.group.eq(g), 'Count'].to_numpy() for g in GROUPS]
                estimable = min(len(sham), len(dbs)) >= 2
                sd = np.sqrt(((len(sham) - 1) * sham.var(ddof=1) + (len(dbs) - 1) * dbs.var(ddof=1)) / (len(sham) + len(dbs) - 2)) if estimable else np.nan
                results.append(dict(target=target, Region=region, n_Sham=len(sham), n_DBS=len(dbs), P=stats.ttest_ind(dbs, sham, equal_var=False).pvalue if estimable else np.nan, CohenD=(dbs.mean() - sham.mean()) / sd if sd > 0 else np.nan, Status='estimable' if estimable else 'fewer_than_two_animals_in_a_group'))
            summary = pd.DataFrame(results)
            summary['FDR'] = np.nan
            valid = summary.P.notna()
            summary.loc[valid, 'FDR'] = stats.false_discovery_control(summary.loc[valid, 'P'].to_numpy())
            summary.to_csv(destination('statistics', f'5{panel_count}_statistics.csv'), index=False)
            fig, ax = plt.subplots(figsize=(8, 3.5))
            for j, group in enumerate(GROUPS):
                for i, region in enumerate(REGIONS):
                    v = frame.loc[frame.group.eq(group) & frame.region.eq(region), 'Count'].to_numpy()
                    pos = i + (j - 0.5) * 0.36
                    ax.bar(pos, v.mean() if len(v) else np.nan, yerr=stats.sem(v) if len(v) > 1 else 0, capsize=2, width=0.33, color=COLORS[j], edgecolor='black', label=group if i == 0 else None)
                    ax.scatter(pos + np.linspace(-0.06, 0.06, len(v)), v, c='white', edgecolors='grey', s=16, zorder=3)
                    q = float(summary.loc[summary.Region.eq(region), 'FDR'].iloc[0])
                    if j == 1 and q < 0.05:
                        ax.text(i, frame.loc[frame.region.eq(region), 'Count'].max() * 1.12, f'q={q:.2g}', ha='center', fontsize=7)
            ax.set(xticks=range(10), xticklabels=REGIONS, ylabel=f'Animal mean NeuN+{target}+ count per section')
            ax.tick_params(axis='x', rotation=45)
            ax.spines[['top', 'right']].set_visible(False)
            ax.legend(frameon=False)
            fig.tight_layout()
            for ext in ['png']:
                fig.savefig(destination('figures', f'5{panel_count}.{ext}'), dpi=400)
            plt.close(fig)
            corr = summary[summary.Region.ne('ATN') & summary.CohenD.notna()].merge(sc2[['Region', 'SC2']], on='Region', validate='one_to_one')
            corr[['Region', 'CohenD', 'SC2', 'n_Sham', 'n_DBS']].to_csv(destination('sourcedata', f'5{panel_corr}.csv'), index=False)
            rho, p = stats.spearmanr(corr.CohenD, corr.SC2)
            all_stats.append(dict(panel=f'5{panel_corr}', rho=rho, P=p, n_regions=len(corr), unit='region; effect from animal mean counts'))
            fig, ax = plt.subplots(figsize=(4, 3.5))
            ax.scatter(corr.CohenD, corr.SC2, color='#65adc0', s=35)
            x = corr.CohenD.to_numpy()
            y = corr.SC2.to_numpy()
            coef = np.polyfit(x, y, 1)
            xx = np.linspace(x.min(), x.max(), 100)
            ax.plot(xx, np.polyval(coef, xx), '--', c='grey')
            for row in corr.itertuples():
                ax.annotate(row.Region, (row.CohenD, row.SC2), xytext=(3, 3), textcoords='offset points', fontsize=7)
            ax.set(xlabel=f"{target} count Cohen's d (DBS - Sham)", ylabel='2nd structural connectivity', title=f'Spearman rho={rho:.2f}, P={p:.3g}')
            ax.spines[['top', 'right']].set_visible(False)
            fig.tight_layout()
            for ext in ['png']:
                fig.savefig(destination('figures', f'5{panel_corr}.{ext}'), dpi=400)
            plt.close(fig)
        pd.DataFrame(all_stats).to_csv(destination('statistics', '5gi_correlations.csv'), index=False)
        print(pd.DataFrame(all_stats).to_string(index=False))
    main()
    return locals()

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
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
    code_base = Path(__file__).resolve().parents[2]
    if OUTPUT.resolve() == code_base or code_base in OUTPUT.resolve().parents:
        raise ValueError('Output must be outside the code package')
    ROOT = DATA / 'results/260814_Rsrc1-Sf3b1 IF/260929_Representative'
    REGIONS = ['ATN', 'CP', 'RSPd', 'RT', 'RE', 'GPe', 'PF', 'DG', 'LSc', 'TRS']

    def destination(kind, name):
        path = OUTPUT / kind / 'extended data figure 14' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def counts(qc_file):
        original = inputs['original_1'][str(ROOT / 'V4/analysis/All_slices_NeuN_Target_Counts_Long_V4.csv')]
        columns = ['target', 'group', 'mouse_id', 'Sample_ID', 'region', 'region_present', 'ROI_area_mm2']
        left = original[columns + ['n_NeuN_Target']].rename(columns={'n_NeuN_Target': 'count'}).assign(hemisphere='Left')
        paths = sorted(inputs['paths_2'])
        if not paths:
            raise FileNotFoundError('No original right-hemisphere count tables found')
        right = pd.concat([inputs['right_3'][str(p)] for p in paths], ignore_index=True)
        data = pd.concat([left, right], ignore_index=True)
        keep = data.region_present.astype(str).str.lower().eq('true') & data.region.isin(REGIONS) & data.group.isin(['DBS', 'Sham'])
        data = data.loc[keep, columns + ['count', 'hemisphere']].copy()
        assert not data.duplicated(['target', 'Sample_ID', 'mouse_id', 'region', 'hemisphere']).any()
        assert data['count'].notna().all() and data['count'].ge(0).all()
        if qc_file:
            exclusions = inputs['exclusions_4'][str(qc_file)]
            keys = ['target', 'Sample_ID', 'region']
            assert set(keys + ['reason']).issubset(exclusions.columns)
            assert not exclusions.duplicated(keys).any()
            assert exclusions.reason.notna().all() and exclusions.reason.astype(str).str.strip().ne('').all()
            merged = exclusions.merge(data[keys].drop_duplicates(), on=keys, how='left', indicator=True)
            assert merged._merge.eq('both').all()
            data = data.merge(exclusions[keys + ['reason']], on=keys, how='left', validate='many_to_one')
            data[data.reason.notna()].to_csv(destination('statistics', 'explicit_qc_exclusions.csv'), index=False)
            data = data[data.reason.isna()].drop(columns='reason')
        return data

    def plot_bar(animal, results, levels, labels, target, panel, paired):
        colors = ['#e9b482', '#b9d7ea'] if paired else ['#91459a', '#ffdc17']
        fig, ax = plt.subplots(figsize=(7, 3.4))
        for i, region in enumerate(REGIONS):
            sub = animal[animal.region.eq(region)]
            for j, condition in enumerate(levels):
                values = sub.loc[sub.condition.eq(condition), 'count'].to_numpy()
                if not len(values):
                    continue
                x = i + (j - 0.5) * 0.34
                ax.bar(x, values.mean(), 0.3, color=colors[j], edgecolor='grey', lw=0.5, label=labels[j] if i == 0 else None)
                if len(values) > 1:
                    ax.errorbar(x, values.mean(), stats.sem(values), fmt='none', color='grey', elinewidth=0.5, capsize=2)
                ax.scatter(x + np.linspace(-0.05, 0.05, len(values)), values, s=15, facecolors='none', edgecolors='grey', lw=0.6)
            q = results.loc[results.Region.eq(region), 'FDR'].iloc[0]
            if q < 0.05:
                ax.text(i, sub['count'].max() * 1.08, f'q={q:.2g}', ha='center', fontsize=7)
        ax.set(xticks=range(10), xticklabels=REGIONS, ylabel=f'Animal mean NeuN+{target}+ count per section')
        ax.tick_params(axis='x', rotation=45)
        ax.legend(frameon=False, ncol=2)
        ax.set_title('Paired hemispheres within DBS animals' if paired else 'Contralateral hemisphere', fontsize=9)
        fig.tight_layout()
        fig.savefig(destination('figures', f'ed14{panel}.png'), dpi=400)
        plt.close(fig)

    def main():
        parser = argparse.ArgumentParser(description='Animal-level bilateral IF statistics; no default sample exclusions.')
        args = inputs['args_5']
        data = counts(args.qc_exclusions)
        imaging = inputs['imaging_6'][str(DATA / 'results/260814_Rsrc1-Sf3b1 IF/260819_Sf3b1/results/tables/imaging_metrics_13_regions.csv')]
        plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
        correlations = []
        for target, paired, bar, scatter in [('Sf3b1', False, 'b', 'c'), ('Rsrc1', False, 'd', 'e'), ('Sf3b1', True, 'f', 'g'), ('Rsrc1', True, 'h', 'i')]:
            frame = data[data.target.eq(target)].copy()
            if paired:
                frame = frame[frame.group.eq('DBS')]
                pivot = frame.pivot(index=['Sample_ID', 'mouse_id', 'region'], columns='hemisphere', values='count')
                pivot = pivot.dropna(subset=['Left', 'Right'])
                source = pivot.reset_index().melt(id_vars=['Sample_ID', 'mouse_id', 'region'], value_vars=['Left', 'Right'], var_name='condition', value_name='count')
                levels = ['Left', 'Right']
                labels = ['DBS_I', 'DBS_C']
                a, b = levels
            else:
                source = frame[frame.hemisphere.eq('Right')].rename(columns={'group': 'condition'})
                levels = ['Sham', 'DBS']
                labels = levels
                a, b = ('DBS', 'Sham')
            source = source[['Sample_ID', 'mouse_id', 'region', 'condition', 'count']].assign(target=target)
            source.to_csv(destination('sourcedata', f'ed14{bar}_sections.csv'), index=False)
            animal = source.groupby(['mouse_id', 'region', 'condition', 'target'], as_index=False).agg(count=('count', 'mean'), n_sections=('Sample_ID', 'size'))
            animal.to_csv(destination('sourcedata', f'ed14{bar}.csv'), index=False)
            records = []
            for region in REGIONS:
                sub = animal[animal.region.eq(region)]
                if paired:
                    pair = sub.pivot(index='mouse_id', columns='condition', values='count').reindex(columns=[a, b]).dropna()
                    va, vb = (pair[a].to_numpy(), pair[b].to_numpy())
                else:
                    va, vb = [sub.loc[sub.condition.eq(c), 'count'].to_numpy() for c in [a, b]]
                p = effect = np.nan
                if min(len(va), len(vb)) >= 2:
                    p = float((stats.ttest_rel(va, vb) if paired else stats.ttest_ind(va, vb, equal_var=False)).pvalue)
                    sd = np.std(va - vb, ddof=1) if paired else np.sqrt(((len(va) - 1) * np.var(va, ddof=1) + (len(vb) - 1) * np.var(vb, ddof=1)) / (len(va) + len(vb) - 2))
                    effect = (va.mean() - vb.mean()) / sd if sd > 0 else np.nan
                records.append(dict(Region=region, P=p, effect=effect, n_a=len(va), n_b=len(vb), Effect_definition='paired Cohen dz (Left-Right)' if paired else 'pooled-SD Cohen d (DBS-Sham)'))
            results = pd.DataFrame(records)
            results['FDR'] = np.nan
            valid = results.P.notna()
            results.loc[valid, 'FDR'] = stats.false_discovery_control(results.loc[valid, 'P'].to_numpy())
            results.to_csv(destination('statistics', f'ed14{bar}_statistics.csv'), index=False)
            plot_bar(animal, results, levels, labels, target, bar, paired)
            corr = results[results.Region.ne('ATN') & results.effect.notna()].merge(imaging[['Region', 'SC2']], on='Region', validate='one_to_one')
            corr[['Region', 'effect', 'SC2', 'n_a', 'n_b', 'Effect_definition']].to_csv(destination('sourcedata', f'ed14{scatter}.csv'), index=False)
            rho, p = stats.spearmanr(corr.effect, corr.SC2)
            correlations.append(dict(panel=f'ed14{scatter}', rho=rho, P=p, n_regions=len(corr)))
            fig, ax = plt.subplots(figsize=(4, 3.6))
            ax.scatter(corr.effect, corr.SC2, c='#7bb3c5', s=30)
            x = corr.effect.to_numpy()
            y = corr.SC2.to_numpy()
            xx = np.linspace(x.min(), x.max(), 100)
            ax.plot(xx, np.polyval(np.polyfit(x, y, 1), xx), '--', color='grey', lw=0.8)
            for row in corr.itertuples():
                ax.annotate(row.Region, (row.effect, row.SC2), xytext=(3, 3), textcoords='offset points', fontsize=7)
            ax.set(xlabel=f'{target}: ' + ('paired Cohen dz (DBS_I - DBS_C)' if paired else 'Cohen d (DBS - Sham)'), ylabel='2nd structural connectivity', title=f'Spearman rho={rho:.2f}; P={p:.3g}')
            fig.tight_layout()
            fig.savefig(destination('figures', f'ed14{scatter}.png'), dpi=400)
            plt.close(fig)
        correlations = pd.DataFrame(correlations)
        correlations['FDR'] = stats.false_discovery_control(correlations.P.to_numpy())
        correlations.to_csv(destination('statistics', 'ed14_correlations.csv'), index=False)
        print(correlations.to_string(index=False))
    main()
    return locals()

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    REGIONS = inputs['REGIONS_1']
    GROUPS = inputs['GROUPS_2']
    COLORS = inputs['COLORS_3']
    group_tests = inputs['group_tests_4']

    def panel_a(context):
        table = inputs['table_5'][str(context.data / 'results/Table/Pathway_Circuit_Integration/1_Pathway_Stats_Master.csv')]
        table = table[table.Pathway.eq('RNA Splicing (GO:0008380)')].copy()
        if table.Region.duplicated().any():
            raise ValueError('Expected one pathway score per region')
        rows = table.melt(id_vars=['Region', 'Pathway'], value_vars=['AUCell_' + g for g in GROUPS], var_name='Group', value_name='AUCell_score')
        rows['Group'] = rows.Group.str.replace('AUCell_', '', regex=False)
        rows['unit'] = 'AUCell score; regional mean'
        context.source_csv(rows, 'a')
        palette = {'CP': '#8d7acc', 'RSPd': '#ff1a71', 'RT': '#ff6600', 'RE': '#1340ff', 'GPe': '#b199ff', 'LH': '#c68105', 'PF': '#0a4093', 'CM': '#08306d', 'DG': '#16bcbc', 'LSc': '#3283fe', 'MH': '#faa307', 'TRS': '#7609b1', 'ATN': '#1460ff'}
        fig, ax = plt.subplots(figsize=(2.7, 3.2))
        labels = []
        for region, frame in rows.groupby('Region', sort=True):
            values = frame.set_index('Group').loc[GROUPS, 'AUCell_score'].to_numpy()
            ax.plot(range(3), values, 'o-', lw=0.7, ms=3, color=palette.get(region, '#777777'))
            labels.append((values[-1], region))
        span = rows.AUCell_score.max() - rows.AUCell_score.min()
        previous = -np.inf
        for value, region in sorted(labels):
            y = max(value, previous + span * 0.055)
            ax.annotate(region, (2, value), xytext=(2.14, y), fontsize=6, color=palette.get(region, '#777777'), arrowprops=dict(arrowstyle='-', color=palette.get(region, '#777777'), lw=0.35))
            previous = y
        ax.set(xticks=range(3), xticklabels=GROUPS, ylabel='RNA splicing AUCell score', xlim=(-0.2, 2.65))
        ax.set_title(f'{table.Region.nunique()} regional means; descriptive', fontsize=7)
        context.save(fig, 'a')

    def panel_b(context):
        cells = context.cells('cell_level_metrics_with_saline.csv')
        cells['spliced_unspliced_ratio'] = cells.spliced_umi / (cells.unspliced_umi + 1.0)
        columns = ['CB', 'Cell_ID', 'Sample', 'donor_id', 'Group', 'final_region', 'spliced_umi', 'unspliced_umi', 'spliced_unspliced_ratio']
        context.source_csv(cells[columns].assign(unit='UMI count; spliced/(unspliced+1) dimensionless'), 'b', '_cells')
        donor = cells.groupby(['donor_id', 'Group', 'final_region']).spliced_unspliced_ratio.agg(['mean', 'count']).reset_index()
        donor = donor.rename(columns={'mean': 'spliced_unspliced_ratio', 'count': 'n_cells'})
        context.source_csv(donor.assign(unit='donor mean of cell ratios'), 'b')
        tests = group_tests(donor, 'spliced_unspliced_ratio', [('Sham', 'Saline'), ('DBS', 'Sham')], context.args.seed)
        tests.to_csv(context.analysis / 'ed11b_donor_tests.csv', index=False)
        fig, ax = plt.subplots(figsize=(8.8, 2.8))
        maximum = donor.spliced_unspliced_ratio.max()
        minimum = donor.spliced_unspliced_ratio.min()
        span = max(maximum - minimum, 0.1)
        for j, group in enumerate(GROUPS):
            for i, region in enumerate(REGIONS):
                values = donor[(donor.Group == group) & (donor.final_region == region)].spliced_unspliced_ratio.to_numpy()
                pos = i + (j - 1) * 0.26
                if not len(values):
                    continue
                ax.bar(pos, values.mean(), 0.23, color=COLORS[group], edgecolor='#555555', lw=0.4, yerr=values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 else None, capsize=1.3, error_kw={'elinewidth': 0.5})
                offsets = np.linspace(-0.055, 0.055, len(values))
                ax.scatter(pos + offsets, values, s=8, facecolors='white', edgecolors='#333333', lw=0.35, zorder=3)
            ax.plot([], [], color=COLORS[group], lw=5, label=group)
        for i, region in enumerate(REGIONS):
            for j, comparison in enumerate(['Sham vs Saline', 'DBS vs Sham']):
                q = tests.loc[(tests.region == region) & (tests.comparison == comparison), 'FDR'].iloc[0]
                text = 'NA' if not np.isfinite(q) else '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'n.s.'
                left, right = (i - 0.26, i) if j == 0 else (i, i + 0.26)
                y = maximum + span * (0.13 + 0.15 * j)
                ax.plot([left, right], [y, y], color='#444444', lw=0.45)
                ax.text((left + right) / 2, y + 0.01 * span, text, ha='center', fontsize=5)
        ax.set(xticks=range(len(REGIONS)), xticklabels=REGIONS, ylabel='Spliced / (Unspliced + 1) UMI', ylim=(0, maximum + span * 0.49))
        ax.legend(frameon=False, ncol=3, loc='upper right', bbox_to_anchor=(1, 1.23))
        ax.set_title('Animal means ± SEM; points are donors; BH-adjusted tests', fontsize=7, loc='left')
        fig.tight_layout()
        context.save(fig, 'b')

    def panel_c(context):
        if context.args.donor_psi:
            psi = inputs['psi_6'][str(context.input(context.args.donor_psi))]
            required = {'event_id', 'donor_id', 'group', 'region', 'PSI'}
            if not required <= set(psi) or psi.duplicated(['donor_id', 'region', 'event_id']).any():
                raise ValueError('Expected unique original donor event PSI')
            if 'read_origin' not in psi or not psi.read_origin.eq('original_STAR_unique').all():
                raise ValueError('Formal panel c requires provenance-tagged original donor PSI')
            records = []
            for region, frame in psi.groupby('region'):
                matrix = frame.pivot(index=['donor_id', 'group'], columns='event_id', values='PSI')
                sham = matrix.loc[matrix.index.get_level_values('group') == 'Sham_I']
                for (donor, group), values in matrix.iterrows():
                    if group not in ['Sham_I', 'DBS_I']:
                        continue
                    reference = sham.loc[sham.index.get_level_values('donor_id') != donor].mean(axis=0)
                    valid = values.notna() & reference.notna()
                    records.append(dict(donor_id=donor, Group=group.replace('_I', ''), final_region=region, M_remodel=(values[valid] - reference[valid]).abs().mean(), n_events=int(valid.sum())))
            donor = pd.DataFrame(records)
            title = 'Original donor PSI; leave-one-donor-out Sham reference'
            error = 'SD'
            context.source_csv(psi, 'c', '_PSI')
        elif context.args.legacy_c_exploratory:
            cells = context.cells('cell_level_5_splicing_metrics.csv')
            cells = cells[cells.Group.isin(['Sham', 'DBS'])]
            context.source_csv(cells[['CB', 'Cell_ID', 'Sample', 'donor_id', 'Group', 'final_region', 'M_remodel']].assign(input_status='legacy_neighbor_aggregated_PSI; exploratory'), 'c', '_cells')
            donor = cells.groupby(['donor_id', 'Group', 'final_region']).M_remodel.agg(['mean', 'count']).reset_index().rename(columns={'mean': 'M_remodel', 'count': 'n_cells'})
            title = 'Legacy aggregated PSI; exploratory donor summaries'
            error = 'SD'
        else:
            raise FileNotFoundError('Panel c requires DBSEQ_DONOR_PSI_FILE from original junction reads. Use --legacy-c-exploratory only for the labeled descriptive historical panel.')
        context.source_csv(donor.assign(unit='mean absolute PSI difference', interpretation=title), 'c')
        fig, ax = plt.subplots(figsize=(7.3, 2.9))
        for j, group in enumerate(['Sham', 'DBS']):
            for i, region in enumerate(REGIONS):
                values = donor[(donor.Group == group) & (donor.final_region == region)].M_remodel.dropna().to_numpy()
                if not len(values):
                    continue
                pos = i + (j - 0.5) * 0.34
                ax.bar(pos, values.mean(), 0.3, color=COLORS[group], edgecolor='#555555', lw=0.4, yerr=values.std(ddof=1) if len(values) > 1 else None, capsize=1.5, error_kw={'elinewidth': 0.5})
                ax.scatter(pos + np.linspace(-0.04, 0.04, len(values)), values, s=8, facecolors='white', edgecolors='#333333', lw=0.35, zorder=3)
            ax.plot([], [], color=COLORS[group], lw=5, label=group)
        ax.set(xticks=range(len(REGIONS)), xticklabels=REGIONS, ylabel='Splicing remodeling magnitude', ylim=(0, None))
        ax.set_title(title + '\nMean ± ' + error + '; points are donor summaries', fontsize=7, loc='left')
        ax.legend(frameon=False, ncol=2)
        fig.tight_layout()
        context.save(fig, 'c')
    return locals()

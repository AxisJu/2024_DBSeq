def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import textwrap
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy import stats
    REGIONS = inputs['REGIONS_1']
    DIVISION = inputs['DIVISION_2']
    bh = inputs['bh_3']

    def event_tables(context):
        raw = inputs['raw_4']
        sig = raw[(raw.FDR < 0.05) & (raw.delta_PSI.abs() >= context.args.delta_psi)].copy()
        common = sig.groupby(['event_id', 'gene_name', 'event_type']).agg(n_regions=('region', 'nunique'), min_fdr=('FDR', 'min'), regions=('region', lambda x: ';'.join(sorted(set(x))))).reset_index()
        return (raw, sig, common)

    def panel_d(context):
        raw, sig, common = event_tables(context)
        counts = common.n_regions.value_counts().reindex(range(1, len(REGIONS) + 1), fill_value=0)
        rows = sig.merge(common[['event_id', 'n_regions']], on='event_id', how='left', validate='many_to_one')
        context.source_csv(rows, 'd')
        fig, ax = plt.subplots(figsize=(4.4, 3.3))
        ax.bar(counts.index, counts.values, color=['#b0bec5'] + ['#37474f'] * (len(REGIONS) - 1))
        ax.set(xlabel='Number of regions sharing a DAS event', ylabel='DAS event count', xticks=range(1, len(REGIONS) + 1, 2))
        if len(common):
            specific = common[common.regions.eq('ATN')].sort_values(['min_fdr', 'event_id']).drop_duplicates('gene_name').head(5)
            shared = common[common.n_regions > 1].sort_values(['n_regions', 'min_fdr', 'event_id'], ascending=[False, True, True]).drop_duplicates('gene_name').head(5)
            labels = ['Top ATN-specific'] + [str(r.gene_name) + ' (' + r.event_type + ')' for r in specific.itertuples()]
            labels += ['', 'Top shared'] + [str(r.gene_name) + f' ({r.n_regions}/{len(REGIONS)})' for r in shared.itertuples()]
            ax.text(0.43, 0.98, '\n'.join(labels), transform=ax.transAxes, va='top', fontsize=6)
        else:
            ax.text(0.5, 0.5, 'No events meet the donor-level criteria', ha='center', transform=ax.transAxes, fontsize=7)
        ax.set_title(f'Donor permutation BH < 0.05; |ΔPSI| ≥ {context.args.delta_psi:g}', fontsize=7)
        fig.tight_layout()
        context.save(fig, 'd')

    def correlation(frame, metric, partial):
        x = frame.DAS_count.to_numpy(float)
        y = frame[metric].to_numpy(float)
        if len(x) < 4 or np.ptp(x) == 0 or np.ptp(y) == 0:
            return (np.nan, np.nan)
        if partial:
            rank_x, rank_y, rank_z = map(stats.rankdata, [x, y, frame.ED.to_numpy(float)])
            design = np.column_stack([np.ones(len(x)), rank_z])
            rx = rank_x - design @ np.linalg.lstsq(design, rank_x, rcond=None)[0]
            ry = rank_y - design @ np.linalg.lstsq(design, rank_y, rcond=None)[0]
            if np.linalg.norm(rx) < 1e-12 or np.linalg.norm(ry) < 1e-12:
                return (np.nan, np.nan)
            rho = stats.pearsonr(rx, ry).statistic
            pv = 2 * stats.t.sf(abs(rho) * np.sqrt((len(x) - 3) / max(1 - rho * rho, 1e-15)), len(x) - 3)
        else:
            rho, pv = stats.spearmanr(x, y)
        return (rho, pv)

    def draw_scatter(ax, frame, metric, record):
        x = frame.DAS_count.to_numpy(float)
        y = frame[metric].to_numpy(float)
        if len(x) >= 3 and np.ptp(x) > 0:
            fit = stats.linregress(x, y)
            xx = np.linspace(x.min(), x.max(), 100)
            ax.plot(xx, fit.intercept + fit.slope * xx, '--', color='#444444', lw=0.7)
        for row in frame.itertuples():
            value = getattr(row, metric)
            ax.scatter(row.DAS_count, value, color=DIVISION.get(row.Division, '#777777'), s=24, edgecolor='white', lw=0.3)
            ax.annotate(row.Region, (row.DAS_count, value), xytext=(3, 3), textcoords='offset points', fontsize=6)
        kind = 'Partial ρ (ED adjusted)' if record['partial'] else 'Spearman ρ'
        label = f"{kind} = {record['rho']:.2f}\nBH q = {record['FDR']:.3g}" if np.isfinite(record['rho']) else 'Correlation not estimable'
        ax.text(0.03, 0.98, label, transform=ax.transAxes, va='top', fontsize=6)
        ax.set(xlabel='Donor-level DAS event count', ylabel=metric)

    def panels_ef(context):
        raw, sig, common = event_tables(context)
        frame = inputs['frame_5'][str(context.data / 'results/Table/imaging_metrics_8metrics_12_regions.csv')]
        frame = frame[frame.Region.ne('ATN') & frame.Region.isin(REGIONS)].copy()
        if frame.Region.duplicated().any():
            raise ValueError('Imaging table has duplicate regions')
        available = set(raw.loc[raw.P.notna(), 'region'])
        missing = set(frame.Region) - available
        if missing:
            raise ValueError('No donor-level DAS tests for imaging regions: ' + ', '.join(sorted(missing)))
        count = sig.groupby('region').event_id.nunique()
        frame['DAS_count'] = frame.Region.map(count).fillna(0).astype(int)
        rows = []
        for metric in ['SC2', 'SPE', 'DE', 'NE', 'CMY']:
            partial = metric != 'SC2'
            rho, pv = correlation(frame.dropna(subset=['DAS_count', metric, 'ED']), metric, partial)
            rows.append(dict(metric=metric, partial=partial, rho=rho, P=pv, n_regions=len(frame)))
        tests = pd.DataFrame(rows)
        tests['FDR'] = bh(tests.P)
        tests.to_csv(context.analysis / 'ed11ef_regional_correlations.csv', index=False)
        if 'e' in context.args.panels:
            context.source_csv(frame[['Region', 'Division', 'DAS_count', 'SC2']].assign(unit='event count; SC2'), 'e')
            fig, ax = plt.subplots(figsize=(3.2, 3.3))
            draw_scatter(ax, frame, 'SC2', tests.iloc[0])
            fig.tight_layout()
            context.save(fig, 'e')
        if 'f' in context.args.panels:
            context.source_csv(frame[['Region', 'Division', 'DAS_count', 'ED', 'SPE', 'DE', 'NE', 'CMY']].assign(unit='event count; ED micrometre; other network metrics dimensionless'), 'f')
            fig, axes = plt.subplots(2, 2, figsize=(6.4, 5.8))
            for ax, metric in zip(axes.flat, ['SPE', 'DE', 'NE', 'CMY']):
                draw_scatter(ax, frame, metric, tests[tests.metric == metric].iloc[0])
            fig.tight_layout()
            context.save(fig, 'f')

    def panel_g(context):
        raw, sig, common = event_tables(context)
        terms = {}
        handle = inputs['handle_6'][str(context.input(context.args.go_gmt))]
        for line in handle:
            term, description, *genes = line.rstrip().split('\t')
            terms[term] = set(genes)
        rows = []
        background_rows = []
        for region, frame in raw.groupby('region'):
            universe = set(frame.loc[frame.P.notna(), 'gene_name'].dropna().astype(str))
            universe -= {'', 'Unknown', 'nan'}
            query = set(sig.loc[sig.region == region, 'gene_name'].dropna().astype(str)) & universe
            background_rows.extend((dict(region=region, gene=gene, is_DAS_gene=gene in query) for gene in sorted(universe)))
            if not universe:
                continue
            for term, members in terms.items():
                eligible = members & universe
                if not eligible:
                    continue
                hit = eligible & query
                pv = stats.hypergeom.sf(len(hit) - 1, len(universe), len(eligible), len(query)) if query else 1.0
                rows.append(dict(region=region, term=term, overlap_count=len(hit), query_size=len(query), term_size=len(eligible), universe_size=len(universe), P=pv, overlap_genes=';'.join(sorted(hit)), Jaccard=len(hit) / len(eligible | query) if eligible | query else 0.0))
        table = pd.DataFrame(rows)
        if table.empty:
            raise ValueError('No tested genes overlap the mouse GO BP background')
        table['FDR'] = bh(table.P)
        table.to_csv(context.analysis / 'ed11g_all_region_term_tests.csv', index=False)
        context.source_csv(pd.DataFrame(background_rows), 'g', '_tested_genes')
        candidates = table[table.overlap_count > 0].sort_values(['FDR', 'P', 'term']).groupby('region', sort=False).head(2)
        display = candidates.groupby('term').FDR.min().sort_values().head(15).index.tolist()
        show = table[table.term.isin(display)].copy()
        context.source_csv(show, 'g')
        fig, ax = plt.subplots(figsize=(8.5, max(3.2, len(display) * 0.26 + 1.2)))
        if len(display):
            positions = {term: i for i, term in enumerate(display)}
            points = ax.scatter(show.region.map({r: i for i, r in enumerate(REGIONS)}), show.term.map(positions), s=show.overlap_count * 10 + 6, c=show.Jaccard, cmap='YlGnBu', edgecolors='#37474f', lw=0.4)
            for row in show.itertuples():
                if row.FDR < 0.05:
                    ax.text(REGIONS.index(row.region), positions[row.term], '*', ha='center', va='center', fontsize=7)
            labels = [textwrap.fill(term.replace('GOBP_', '').replace('_', ' ').lower(), 46) for term in display]
            ax.set(yticks=range(len(display)), yticklabels=labels)
            ax.invert_yaxis()
            fig.colorbar(points, ax=ax, label='Jaccard', fraction=0.025, pad=0.025)
        else:
            ax.text(0.5, 0.5, 'No DAS genes overlap the tested GO BP terms', ha='center', transform=ax.transAxes)
        ax.set(xticks=range(len(REGIONS)), xticklabels=REGIONS)
        ax.set_title('DAS gene enrichment; region-specific tested-gene background; * BH < 0.05', fontsize=7)
        fig.tight_layout()
        context.save(fig, 'g')
    return locals()

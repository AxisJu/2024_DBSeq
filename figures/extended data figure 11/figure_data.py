def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    from pathlib import Path
    import os
    import subprocess
    import sys
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy import stats
    REGIONS = ['ATN', 'CP', 'RSPd', 'RT', 'RE', 'GPe', 'LH', 'PF', 'CM', 'DG', 'LSc', 'MH', 'TRS']
    GROUPS = ['Saline', 'Sham', 'DBS']
    COLORS = dict(Saline='#18a799', Sham='#ac00cc', DBS='#ffd700')
    DIVISION = {'CNU': '#9bd5f4', 'Isocortex': '#18a799', 'TH': '#e0abce', 'HY': '#e74538', 'HPF': '#84c551'}
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': 0.5, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5})

    def bh(values):
        values = np.asarray(values, dtype=float)
        result = np.full(values.shape, np.nan)
        valid = np.flatnonzero(np.isfinite(values))
        order = valid[np.argsort(values[valid])]
        if len(order):
            result[order] = np.minimum(1, np.minimum.accumulate((values[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1])
        return result

    def group_tests(donors, value, comparisons, seed):
        rows = []
        rng = np.random.default_rng(seed)
        for region in REGIONS:
            frame = donors[donors.final_region == region]
            for b, a in comparisons:
                x = frame.loc[frame.Group == b, value].dropna().to_numpy()
                y = frame.loc[frame.Group == a, value].dropna().to_numpy()
                p = np.nan
                if min(len(x), len(y)) >= 2:
                    p = stats.permutation_test((x, y), lambda u, v: np.mean(u) - np.mean(v), vectorized=False, permutation_type='independent', alternative='two-sided', n_resamples=9999, random_state=rng).pvalue
                rows.append(dict(region=region, comparison=b + ' vs ' + a, n_b=len(x), n_a=len(y), P=p))
        result = pd.DataFrame(rows)
        result['FDR'] = bh(result.P)
        return result

    class Context:

        def __init__(self, args):
            self.args = args
            if not os.environ.get('DBSEQ_DATA_ROOT') or not os.environ.get('DBSEQ_OUTPUT_ROOT'):
                raise ValueError('Set DBSEQ_DATA_ROOT and DBSEQ_OUTPUT_ROOT')
            self.data = Path(os.environ['DBSEQ_DATA_ROOT']).resolve()
            self.output = Path(os.environ['DBSEQ_OUTPUT_ROOT']).resolve()
            self.package = Path(__file__).resolve().parents[2]
            if self.output == self.package or self.package in self.output.parents:
                raise ValueError('DBSEQ_OUTPUT_ROOT must be outside the publication code package')
            self.fig = self.output / 'figures/extended data figure 11'
            self.source = self.output / 'sourcedata/extended data figure 11'
            self.analysis = self.output / 'analysis_checks/extended data figure 11'
            for path in (self.fig, self.source, self.analysis):
                path.mkdir(parents=True, exist_ok=True)
            self.base = self.data / 'results/260813 DOLPHIN'
            self.das = None

        def input(self, value):
            path = Path(value).expanduser()
            return path.resolve() if path.is_absolute() else self.data / path

        def save(self, fig, panel):
            fig.savefig(self.fig / f'ed11{panel}.png', dpi=400, bbox_inches='tight')
            plt.close(fig)

        def source_csv(self, frame, panel, suffix=''):
            frame.to_csv(self.source / f'ed11{panel}{suffix}.csv', index=False)

        def cells(self, filename):
            if self.args.donor_metadata:
                path = self.input(self.args.donor_metadata)
            else:
                path = self.output / '600_Others/260909-DOLPHIN/metadata_donor.csv'
                entry = self.package / 'analysis/600_Others/260909-DOLPHIN/00_metadata.py'
                inputs['data_1']
            metadata = inputs['metadata_2'][str(path)]
            keys = ['CB', 'Sample']
            if metadata.duplicated(keys).any():
                raise ValueError('Duplicate donor mapping keys')
            cells = inputs['cells_3'][str(self.base / filename)]
            if cells.duplicated(keys).any():
                raise ValueError('Duplicate cell observation keys')
            cells = cells.merge(metadata[keys + ['Cell_ID', 'donor_id']], on=keys, how='left', validate='one_to_one')
            if cells.donor_id.isna().any():
                raise ValueError('Some observations have no biological donor mapping')
            cells = cells[cells.final_region.isin(REGIONS) & cells.Group.isin(GROUPS)].copy()
            if (cells.groupby('donor_id').Group.nunique() > 1).any():
                raise ValueError('A donor occurs in multiple conditions; this design is not paired')
            return cells

        def load_das(self):
            if self.das is not None:
                return self.das
            if not self.args.das_file:
                raise FileNotFoundError('Panels d-g require DBSEQ_DAS_FILE or --das-file from 07_donor_das.py. Supply donor-permutation results from original junction counts. Original unaggregated STAR junctions/BAMs are required; use --panels ab for available formal panels.')
            frame = inputs['frame_4'][str(self.input(self.args.das_file))]
            required = {'event_id', 'gene_name', 'event_type', 'region', 'PSI_Sham', 'PSI_DBS', 'delta_PSI', 'P', 'FDR', 'n_donors_Sham', 'n_donors_DBS', 'test', 'read_origin'}
            if not required <= set(frame) or frame.duplicated(['event_id', 'region']).any():
                raise ValueError('Incomplete or duplicated donor DAS table')
            if not frame.test.eq('donor_permutation').all() or not frame.read_origin.eq('original_STAR_unique').all():
                raise ValueError('Expected donor permutation results from original STAR unique reads')
            tested = frame.P.notna()
            if (frame.loc[tested, ['n_donors_Sham', 'n_donors_DBS']] < 2).any().any():
                raise ValueError('DAS P values require biological replication in each condition')
            if not np.allclose(frame.FDR.to_numpy(), bh(frame.P), equal_nan=True):
                raise ValueError('FDR must use the complete region-by-event testing family')
            frame = frame[frame.region.isin(REGIONS)].copy()
            self.das = frame
            return frame
    return locals()

# Project adaptations of DOLPHIN: https://github.com/mcgilldinglab/DOLPHIN
# MIT License. Copyright (c) 2024 McGill Ding Lab.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import os
    import numpy as np
    import pandas as pd
    import anndata
    from scipy import sparse
    from tqdm import tqdm
    import warnings
    from anndata._core.views import ImplicitModificationWarning
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=ImplicitModificationWarning)

    def _parse_gene_order(val):
        s = str(val)
        try:
            return int(s[4:])
        except Exception:
            digits = ''.join((c for c in s if c.isdigit()))
            return int(digits) if digits else 0

    def _parse_exon_order(val):
        try:
            return int(str(val).split('-')[-1])
        except Exception:
            return 0

    def run_adjacency_matrix_final(out_name: str, out_directory: str='./', batch_size: int=1000):
        print('=' * 65)
        print('>>> [Step 6] Generating Final Cleaned Adjacency Matrix...')
        print('=' * 65)
        final_out_dir = os.path.join(out_directory, 'data')
        adj_file = os.path.join(final_out_dir, f'AdjacencyComp_{out_name}.h5ad')
        fea_file = os.path.join(final_out_dir, f'Feature_{out_name}.h5ad')
        print('Loading Adjacency and Feature matrices into memory (Sparse Mode)...')
        adata_adj = inputs['adata_adj_1'][str(adj_file)]
        if not sparse.isspmatrix_csr(adata_adj.X):
            adata_adj.X = sparse.csr_matrix(adata_adj.X)
        adata_fea_gtf = inputs['adata_fea_gtf_2'][str(fea_file)]
        print('Identifying expressed exons across all cells...')
        fea_var = pd.DataFrame(index=adata_fea_gtf.var.index)
        fea_var['gene_id'] = adata_fea_gtf.var['gene_id'].values
        fea_var['gtf_index'] = fea_var.index
        fea_var['gene_order'] = fea_var['gene_id'].apply(_parse_gene_order)
        fea_var['exon_order'] = fea_var['gtf_index'].apply(_parse_exon_order)
        if sparse.issparse(adata_fea_gtf.X):
            fea_csr = adata_fea_gtf.X.tocsr()
            fea_has_expr = np.zeros(adata_fea_gtf.shape[1], dtype=bool)
            if len(fea_csr.data) > 0:
                fea_has_expr[fea_csr.indices[fea_csr.data > 0]] = True
        else:
            fea_has_expr = (adata_fea_gtf.X > 0).any(axis=0)
        fea_var['has_expr'] = fea_has_expr
        fea_var_order = fea_var.sort_values(by=['gene_order', 'exon_order']).reset_index(drop=True)
        fea_var_order['new_gtf_exon_index'] = fea_var_order.groupby('gene_id', observed=False).cumcount() + 1
        expressed = fea_var_order[fea_var_order['has_expr']]
        exon_keep_dict = expressed.groupby('gene_id', observed=False)['new_gtf_exon_index'].apply(lambda s: set(s.astype(str))).to_dict()
        print('Computing column-wise maximum connectivity across 4.87M+ features...')
        col_max = np.zeros(adata_adj.shape[1], dtype=adata_adj.X.dtype)
        if len(adata_adj.X.data) > 0:
            np.maximum.at(col_max, adata_adj.X.indices, adata_adj.X.data)
        print('Grouping adjacency feature coordinates by gene...')
        var_df = adata_adj.var.copy()
        if 'gene_name' not in var_df.columns:
            var_df['gene_name'] = [str(x)[:str(x).rfind('-')] if '-' in str(x) else str(x) for x in var_df.index]
        if 'gene_id' not in var_df.columns:
            gene_name_to_id = dict(zip(var_df['gene_name'], var_df['gene_id'])) if 'gene_id' in var_df.columns else {}
            var_df['gene_id'] = var_df['gene_name'].map(gene_name_to_id).fillna(var_df['gene_name'])
        gene_ids = var_df['gene_id'].values
        gene_names = var_df['gene_name'].values
        gene_df = pd.DataFrame({'gene_id': gene_ids, 'gene_name': gene_names, 'col_idx': np.arange(len(gene_ids))})
        print('Filtering unexpressed and disconnected exon-pair edges...')
        kept_cols_list = []
        new_var_names_list = []
        new_gene_ids_list = []
        new_gene_names_list = []
        grouped = gene_df.groupby('gene_id', sort=False, observed=False)
        for gene_id, grp in tqdm(grouped, desc='Processing Genes'):
            cols = grp['col_idx'].values
            gene_name = grp['gene_name'].iloc[0]
            adj_size = len(cols)
            n = int(np.sqrt(adj_size))
            if n * n != adj_size or n == 0:
                kept_cols_list.append(cols)
                for idx in range(1, len(cols) + 1):
                    new_var_names_list.append(f'{gene_id}-{idx}')
                    new_gene_ids_list.append(gene_id)
                    new_gene_names_list.append(gene_name)
                continue
            matrix = col_max[cols].reshape((n, n))
            keep_ids = exon_keep_dict.get(gene_id, set())
            removable = np.zeros(n, dtype=bool)
            for i in range(n):
                if np.all(matrix[i, :] == 0) and np.all(matrix[:, i] == 0) and (str(i + 1) not in keep_ids):
                    removable[i] = True
            kept_exons = np.where(~removable)[0]
            m = len(kept_exons)
            if m > 0:
                r_grid, c_grid = np.meshgrid(kept_exons, kept_exons, indexing='ij')
                kept_rel_indices = (r_grid * n + c_grid).flatten()
                kept_cols_list.append(cols[kept_rel_indices])
                for idx in range(1, m * m + 1):
                    new_var_names_list.append(f'{gene_id}-{idx}')
                    new_gene_ids_list.append(gene_id)
                    new_gene_names_list.append(gene_name)
        all_kept_cols = np.concatenate(kept_cols_list)
        print(f'Original features: {adata_adj.shape[1]} -> Cleaned features: {len(all_kept_cols)}')
        print('Slicing sparse matrix and assembling final AnnData...')
        final_X = adata_adj.X[:, all_kept_cols]
        final_var = pd.DataFrame(index=new_var_names_list)
        final_var['gene_id'] = new_gene_ids_list
        final_var['gene_name'] = new_gene_names_list
        final_adata = anndata.AnnData(X=final_X, obs=adata_adj.obs.copy(), var=final_var)
        final_output_path = os.path.join(final_out_dir, f'AdjacencyCompRe_{out_name}.h5ad')
        print(f'Writing final cleaned AnnData to: {final_output_path}')
        final_adata.write_h5ad(final_output_path, compression='gzip')
        print(f'Step 6 Finished Successfully! Final Shape: {final_adata.shape}')
    return locals()

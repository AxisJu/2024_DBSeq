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
    warnings.filterwarnings('ignore')

    def run_adjacency_hvg(out_name: str, out_directory: str='./'):
        print('=' * 65)
        print('>>> [Step 9] Within-Gene Graph Edge Normalization (ALL Genes Retained)...')
        print('=' * 65)
        final_out_dir = os.path.join(out_directory, 'data')
        adj_file = os.path.join(final_out_dir, f'AdjacencyCompRe_{out_name}.h5ad')
        print(f'Loading Cleaned Adjacency matrix: {adj_file} ...')
        adj_anndata = inputs['adj_anndata_1'][str(adj_file)]
        if not sparse.isspmatrix_csr(adj_anndata.X):
            adj_anndata.X = sparse.csr_matrix(adj_anndata.X)
        raw_hvg_file = os.path.join(final_out_dir, f'AdjacencyCompReHvg_{out_name}.h5ad')
        print(f'Writing unnormalized HVG-compatible matrix to: {raw_hvg_file} ...')
        adj_anndata.write_h5ad(raw_hvg_file, compression='gzip')
        print('Pre-computing column-to-gene mapping index...')
        var_df = adj_anndata.var
        if 'gene_id' in var_df.columns:
            gene_series = var_df['gene_id'].astype(str)
        else:
            gene_series = pd.Series([str(x).rsplit('-', 1)[0] for x in var_df.index])
        unique_genes, gene_inverse = np.unique(gene_series.values, return_inverse=True)
        num_unique_genes = len(unique_genes)
        print('Performing vectorized within-gene edge normalization across all cells...')
        X_csr = adj_anndata.X.copy()
        indptr = X_csr.indptr
        indices = X_csr.indices
        data = X_csr.data.copy()
        for i in tqdm(range(X_csr.shape[0]), desc='Normalizing Edges'):
            start, end = (indptr[i], indptr[i + 1])
            if start == end:
                continue
            row_cols = indices[start:end]
            row_data = data[start:end]
            row_genes = gene_inverse[row_cols]
            gene_edge_sums = np.bincount(row_genes, weights=row_data, minlength=num_unique_genes)
            divisors = gene_edge_sums[row_genes]
            mask = divisors > 0
            data[start:end][mask] = row_data[mask] / divisors[mask]
        X_csr.data = data
        new_adata = anndata.AnnData(X=X_csr, obs=adj_anndata.obs.copy(), var=adj_anndata.var.copy())
        edge_file = os.path.join(final_out_dir, f'AdjacencyCompReHvgEdge_{out_name}.h5ad')
        print(f'Writing normalized Edge matrix to: {edge_file} ...')
        new_adata.write_h5ad(edge_file, compression='gzip')
        print(f'Step 9 Finished! Retained {len(unique_genes)} genes | Final Shape: {new_adata.shape}')
    return locals()

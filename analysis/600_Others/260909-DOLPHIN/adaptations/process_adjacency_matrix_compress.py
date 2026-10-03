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
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    import concurrent.futures
    import numpy as np
    import pandas as pd
    import anndata
    from scipy import sparse
    from tqdm import tqdm
    import warnings
    from anndata._core.views import ImplicitModificationWarning
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=ImplicitModificationWarning)
    worker_adj_X = None
    worker_fea_X = None
    worker_gene_fea_cols = None
    worker_gene_adj_cols = None
    worker_gene_sizes = None
    worker_fea_col_to_gene_idx = None
    worker_all_valid_adj_cols = None
    worker_var_df = None
    worker_out_path_dir = None
    worker_obs_name_to_idx = None
    worker_num_genes = 0

    def worker_init(adj_X, fea_X, g_fea_cols, g_adj_cols, g_sizes, fea_col_map, all_valid_adj, v_df, out_dir, obs_to_idx, n_genes):
        global worker_adj_X, worker_fea_X, worker_gene_fea_cols, worker_gene_adj_cols
        global worker_gene_sizes, worker_fea_col_to_gene_idx, worker_all_valid_adj_cols
        global worker_var_df, worker_out_path_dir, worker_obs_name_to_idx, worker_num_genes
        worker_adj_X = adj_X
        worker_fea_X = fea_X
        worker_gene_fea_cols = g_fea_cols
        worker_gene_adj_cols = g_adj_cols
        worker_gene_sizes = g_sizes
        worker_fea_col_to_gene_idx = fea_col_map
        worker_all_valid_adj_cols = all_valid_adj
        worker_var_df = v_df
        worker_out_path_dir = out_dir
        worker_obs_name_to_idx = obs_to_idx
        worker_num_genes = n_genes

    def worker_task(sample_name):
        out_file = os.path.join(worker_out_path_dir, f'{sample_name}.h5ad')
        if inputs['data_1'] and os.path.getsize(out_file) > 1024:
            return True
        try:
            row_idx = worker_obs_name_to_idx[sample_name]
            fea_row = worker_fea_X[row_idx]
            cell_adj_dense = worker_adj_X[row_idx, :].toarray().flatten()
            orig_adj_copy = cell_adj_dense.copy()
            cell_adj_dense[worker_all_valid_adj_cols] = 0
            valid_mask = fea_row.data > 0
            active_fea_indices = fea_row.indices[valid_mask]
            if len(active_fea_indices) > 0:
                gene_mapped = worker_fea_col_to_gene_idx[active_fea_indices]
                valid_mapped = gene_mapped[gene_mapped >= 0]
                if len(valid_mapped) > 0:
                    counts = np.bincount(valid_mapped, minlength=worker_num_genes)
                    candidate_genes = np.nonzero(counts >= 2)[0]
                    cell_fea_dense = fea_row.toarray().flatten()
                    for g_idx in candidate_genes:
                        f_cols = worker_gene_fea_cols[g_idx]
                        a_cols = worker_gene_adj_cols[g_idx]
                        size = worker_gene_sizes[g_idx]
                        temp_fea = cell_fea_dense[f_cols]
                        non_zero = np.nonzero(temp_fea > 0)[0]
                        if len(non_zero) > 1:
                            temp_adj = orig_adj_copy[a_cols]
                            old_adj = temp_adj.reshape((size, size))
                            new_adj = np.zeros((size, size), dtype=old_adj.dtype)
                            r, c = np.triu_indices(len(non_zero), k=1)
                            orig_r = non_zero[r]
                            orig_c = non_zero[c]
                            new_adj[orig_r, orig_c] = old_adj[orig_r, orig_c]
                            adj_rows = non_zero[:-1]
                            adj_cols = non_zero[1:]
                            mask = new_adj[adj_rows, adj_cols] == 0
                            new_adj[adj_rows[mask], adj_cols[mask]] = 1
                            cell_adj_dense[a_cols] = new_adj.flatten()
            new_X = sparse.csr_matrix(cell_adj_dense)
            cell_adata = anndata.AnnData(X=new_X, obs=pd.DataFrame(index=[sample_name]), var=worker_var_df)
            cell_adata.write_h5ad(out_file)
            return True
        except Exception as e:
            print(f'Error on {sample_name}: {e}')
            return False

    def run_adjacency_compression(metadata_path, out_name, out_directory, num_processes=40):
        print('>>> [1/4] Loading Sparse Matrices into RAM...')
        adata_adj_orig = inputs['adata_adj_orig_2'][str(os.path.join(out_directory, 'data', f'Adjacency_{out_name}.h5ad'))]
        if not sparse.isspmatrix_csr(adata_adj_orig.X):
            adata_adj_orig.X = sparse.csr_matrix(adata_adj_orig.X)
        adata_fea_orig = inputs['adata_fea_orig_3'][str(os.path.join(out_directory, 'data', f'Feature_{out_name}.h5ad'))]
        if not sparse.isspmatrix_csr(adata_fea_orig.X):
            adata_fea_orig.X = sparse.csr_matrix(adata_fea_orig.X)
        print('>>> [2/4] Pre-computing Compact Lookup Arrays...')
        fea_gene_ids = adata_fea_orig.var['gene_id'].values
        adj_gene_ids = adata_adj_orig.var['gene_id'].values
        fea_df = pd.DataFrame({'gene_id': fea_gene_ids, 'idx': np.arange(len(fea_gene_ids))})
        adj_df = pd.DataFrame({'gene_id': adj_gene_ids, 'idx': np.arange(len(adj_gene_ids))})
        fea_grouped = fea_df.groupby('gene_id', observed=False)['idx'].apply(np.array).to_dict()
        adj_grouped = adj_df.groupby('gene_id', observed=False)['idx'].apply(np.array).to_dict()
        gene_expr_sum = np.array(adata_fea_orig.X.sum(axis=0)).flatten()
        nonzero_gene_set = {gene for gene, indices in fea_grouped.items() if gene_expr_sum[indices].sum() > 0}
        var_indices = adata_fea_orig.var.index.values
        exon_gene_set = set()
        for gene, indices in fea_grouped.items():
            try:
                if int(var_indices[indices[-1]].split('-')[-1]) > 1:
                    exon_gene_set.add(gene)
            except (ValueError, IndexError):
                pass
        gene_list = sorted(list(nonzero_gene_set & exon_gene_set))
        num_genes = len(gene_list)
        gene_fea_cols = [fea_grouped[g] for g in gene_list]
        gene_adj_cols = [adj_grouped[g] for g in gene_list]
        gene_sizes = np.array([int(np.sqrt(len(cols))) for cols in gene_adj_cols], dtype=np.int32)
        all_valid_adj_cols = np.concatenate(gene_adj_cols)
        fea_col_map = np.full(adata_fea_orig.shape[1], -1, dtype=np.int32)
        for g_idx, cols in enumerate(gene_fea_cols):
            fea_col_map[cols] = g_idx
        print('>>> [3/4] Aligning Metadata and Checking Target Directory...')
        pd_gt = inputs['pd_gt_4'][str(metadata_path)]
        raw_sample_list = list(pd_gt['CB'])
        obs_set = set(adata_adj_orig.obs_names)
        valid_samples = [s for s in raw_sample_list if s in obs_set]
        out_path_dir = os.path.join(out_directory, 'data', 'temp', 'adj_comp_matrix')
        os.makedirs(out_path_dir, exist_ok=True)
        var_df = adata_adj_orig.var.copy()
        obs_name_to_idx = {name: idx for idx, name in enumerate(adata_adj_orig.obs_names)}
        print('=' * 65)
        print(f'Target Valid Multi-Exon Genes : {num_genes}')
        print(f'Matched Cells to Process       : {len(valid_samples)}')
        print(f'Active Worker Processes        : {num_processes}')
        print('=' * 65)
        print('>>> [4/4] Launching High-Throughput Processing Pool...')
        init_args = (adata_adj_orig.X, adata_fea_orig.X, gene_fea_cols, gene_adj_cols, gene_sizes, fea_col_map, all_valid_adj_cols, var_df, out_path_dir, obs_name_to_idx, num_genes)
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes, initializer=worker_init, initargs=init_args) as executor:
            results = list(tqdm(executor.map(worker_task, valid_samples, chunksize=25), total=len(valid_samples), desc='Compressing Cells'))
        print(f'\nExecution Complete! Successfully processed {len(results)} cells.')
    return locals()

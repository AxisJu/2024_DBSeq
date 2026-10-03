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
    import numpy as np
    import pandas as pd
    import anndata
    from scipy import sparse
    import torch
    from torch_geometric.data import Data
    from tqdm import tqdm
    import warnings
    warnings.filterwarnings('ignore')

    def run_model_input(metadata_path: str, out_name: str, out_directory: str='./', gnn_run_num: int=1000, celltypename: str=None):
        print('=' * 65)
        print('>>> [Step 10] Assembling Vectorized PyTorch Geometric Graph Dataset...')
        print('=' * 65)
        final_out_dir = os.path.join(out_directory, 'data')
        save_dir = os.path.join(final_out_dir, 'model_chunks')
        os.makedirs(save_dir, exist_ok=True)
        print('Loading Sparse Matrices into memory...')
        adata_fea = inputs['adata_fea_1'][str(os.path.join(final_out_dir, f'FeatureCompHvg_{out_name}.h5ad'))]
        adata_adj_raw = inputs['adata_adj_raw_2'][str(os.path.join(final_out_dir, f'AdjacencyCompReHvg_{out_name}.h5ad'))]
        adata_adj_edge = inputs['adata_adj_edge_3'][str(os.path.join(final_out_dir, f'AdjacencyCompReHvgEdge_{out_name}.h5ad'))]
        if not sparse.isspmatrix_csr(adata_fea.X):
            adata_fea.X = sparse.csr_matrix(adata_fea.X)
        if not sparse.isspmatrix_csr(adata_adj_raw.X):
            adata_adj_raw.X = sparse.csr_matrix(adata_adj_raw.X)
        if not sparse.isspmatrix_csr(adata_adj_edge.X):
            adata_adj_edge.X = sparse.csr_matrix(adata_adj_edge.X)
        total_cells = adata_fea.shape[0]
        total_nodes = adata_fea.shape[1]
        total_adj_cols = adata_adj_edge.shape[1]
        sample_list = list(adata_fea.obs_names)
        print(f'Total Cells: {total_cells} | Total Exon Nodes: {total_nodes} | Total Edge Features: {total_adj_cols}')
        print('Pre-computing Global Graph Coordinate Mapping Arrays...')
        fea_var = pd.DataFrame(index=adata_fea.var.index)
        fea_var['gene_id'] = adata_fea.var['gene_id'].values
        gene_order = sorted(list(set(fea_var['gene_id'])))
        gene_to_offset = {}
        curr_offset = 0
        fea_gene_counts = fea_var.groupby('gene_id', sort=False, observed=False).size()
        for gid in gene_order:
            gene_to_offset[gid] = curr_offset
            curr_offset += fea_gene_counts.get(gid, 0)
        adj_var = pd.DataFrame(index=adata_adj_edge.var.index)
        adj_var['gene_id'] = adata_adj_edge.var['gene_id'].values
        global_src = np.zeros(total_adj_cols, dtype=np.int64)
        global_dst = np.zeros(total_adj_cols, dtype=np.int64)
        adj_grouped = adj_var.groupby('gene_id', sort=False, observed=False)
        for gid, grp in tqdm(adj_grouped, desc='Mapping Edge Coordinates'):
            cols = grp.index.values if isinstance(grp.index[0], int) else np.arange(len(grp))
            col_indices = np.where(adj_var['gene_id'].values == gid)[0]
            n_cols = len(col_indices)
            size = int(np.sqrt(n_cols))
            if size * size == n_cols and gid in gene_to_offset:
                base_node = gene_to_offset[gid]
                r_rel = np.arange(n_cols) // size
                c_rel = np.arange(n_cols) % size
                global_src[col_indices] = base_node + r_rel
                global_dst[col_indices] = base_node + c_rel
        df_label = inputs['df_label_4'][str(metadata_path)]
        sample_col = 'CB' if 'CB' in df_label.columns else df_label.columns[0]
        df_label.set_index(sample_col, inplace=True)
        if celltypename and celltypename in df_label.columns:
            unique_types = sorted(df_label[celltypename].dropna().unique())
            mapper = {ctype: idx for idx, ctype in enumerate(unique_types)}
            mapped_labels = df_label[celltypename].map(mapper).fillna(0).astype(int)
        else:
            mapped_labels = pd.Series(0, index=df_label.index)
        cell_labels = [mapped_labels.get(cb, 0) for cb in sample_list]
        print(f'Constructing and saving PyG Data chunks (Chunk size = {gnn_run_num})...')
        num_chunks = int(np.ceil(total_cells / gnn_run_num))
        for chunk_idx in range(num_chunks):
            start_i = chunk_idx * gnn_run_num
            end_i = min((chunk_idx + 1) * gnn_run_num, total_cells)
            chunk_data_list = []
            for i in range(start_i, end_i):
                s_name = sample_list[i]
                y_val = cell_labels[i]
                fea_vec = adata_fea.X[i].toarray().flatten()
                cell_x = torch.from_numpy(fea_vec.reshape(-1, 1)).float()
                x_fea_t = torch.from_numpy(fea_vec).float()
                adj_raw_vec = adata_adj_raw.X[i].toarray().flatten()
                x_adj_t = torch.from_numpy(adj_raw_vec).float()
                edge_row = adata_adj_edge.X[i]
                if edge_row.nnz > 0:
                    active_cols = edge_row.indices
                    active_weights = edge_row.data
                    src_nodes = global_src[active_cols]
                    dst_nodes = global_dst[active_cols]
                    edge_index_t = torch.from_numpy(np.stack([src_nodes, dst_nodes], axis=0)).long()
                    edge_attr_t = torch.from_numpy(active_weights.reshape(-1, 1)).float()
                else:
                    edge_index_t = torch.empty((2, 0), dtype=torch.long)
                    edge_attr_t = torch.empty((0, 1), dtype=torch.float32)
                data_obj = Data(x=cell_x, edge_index=edge_index_t, edge_attr=edge_attr_t, y=torch.tensor(y_val, dtype=torch.long), x_fea=x_fea_t, x_adj=x_adj_t, sample_name=s_name)
                chunk_data_list.append(data_obj)
            chunk_path = os.path.join(save_dir, f'model_{out_name}_chunk_{chunk_idx}.pt')
            torch.save(chunk_data_list, chunk_path)
            print(f'[{chunk_idx + 1}/{num_chunks}] Saved chunk ({len(chunk_data_list)} cells) -> {chunk_path}')
        manifest = {'out_name': out_name, 'total_cells': total_cells, 'num_chunks': num_chunks, 'chunk_size': gnn_run_num, 'chunk_dir': save_dir, 'num_nodes': total_nodes, 'num_adj_features': total_adj_cols}
        torch.save(manifest, os.path.join(final_out_dir, f'model_{out_name}_manifest.pt'))
        print('\n' + '=' * 65)
        print('>>> Step 10 Finished Successfully!')
        print(f'Saved {num_chunks} chunks to {save_dir}')
        print('=' * 65)
    return locals()

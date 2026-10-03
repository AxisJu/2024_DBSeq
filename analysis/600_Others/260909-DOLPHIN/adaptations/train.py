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
    os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'
    import pyro
    import torch
    from torch.utils.data import Dataset
    from torch_geometric.data import Data
    from torch_geometric.loader import DataLoader
    import numpy as np
    import pandas as pd
    import anndata
    from scipy import sparse
    from tqdm import tqdm
    import gc
    define_svi = inputs['define_svi_1']
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True

    class DolphinLazyDataset(Dataset):

        def __init__(self, fea_h5ad_path, adj_raw_h5ad_path, adj_edge_h5ad_path, label_series=None):
            print('Loading sparse CSR matrices into memory...')
            self.adata_fea = inputs['self_adata_fea_2'][str(fea_h5ad_path)]
            self.adata_adj_raw = inputs['self_adata_adj_raw_3'][str(adj_raw_h5ad_path)]
            self.adata_adj_edge = inputs['self_adata_adj_edge_4'][str(adj_edge_h5ad_path)]
            self.fea_csr = self.adata_fea.X if sparse.isspmatrix_csr(self.adata_fea.X) else sparse.csr_matrix(self.adata_fea.X)
            self.adj_raw_csr = self.adata_adj_raw.X if sparse.isspmatrix_csr(self.adata_adj_raw.X) else sparse.csr_matrix(self.adata_adj_raw.X)
            self.adj_edge_csr = self.adata_adj_edge.X if sparse.isspmatrix_csr(self.adata_adj_edge.X) else sparse.csr_matrix(self.adata_adj_edge.X)
            self.total_cells = self.fea_csr.shape[0]
            self.num_nodes = self.fea_csr.shape[1]
            self.total_adj_cols = self.adj_edge_csr.shape[1]
            self.labels = label_series.values if label_series is not None else np.zeros(self.total_cells, dtype=np.int64)
            print('Pre-computing Edge Coordinate Maps...')
            fea_var = pd.DataFrame(index=self.adata_fea.var.index)
            fea_var['gene_id'] = self.adata_fea.var['gene_id'].values
            gene_positions = {gid: np.flatnonzero(fea_var['gene_id'].to_numpy() == gid) for gid in fea_var['gene_id'].unique()}
            if not self.adata_fea.obs_names.equals(self.adata_adj_raw.obs_names) or not self.adata_fea.obs_names.equals(self.adata_adj_edge.obs_names):
                raise ValueError('Feature and adjacency cell IDs are not aligned')
            if not self.adata_adj_raw.var_names.equals(self.adata_adj_edge.var_names):
                raise ValueError('Raw and binary adjacency columns are not aligned')
            adj_var = pd.DataFrame(index=self.adata_adj_edge.var.index)
            adj_var['gene_id'] = self.adata_adj_edge.var['gene_id'].values
            self.global_src = np.zeros(self.total_adj_cols, dtype=np.int64)
            self.global_dst = np.zeros(self.total_adj_cols, dtype=np.int64)
            adj_grouped = adj_var.groupby('gene_id', sort=False, observed=False)
            for gid, grp in adj_grouped:
                col_indices = np.where(adj_var['gene_id'].values == gid)[0]
                n_cols = len(col_indices)
                size = int(np.sqrt(n_cols))
                if gid not in gene_positions or size * size != n_cols or len(gene_positions[gid]) != size:
                    raise ValueError('Adjacency block does not match exon features: ' + str(gid))
                positions = gene_positions[gid]
                self.global_src[col_indices] = positions[np.arange(n_cols) // size]
                self.global_dst[col_indices] = positions[np.arange(n_cols) % size]

        def __len__(self):
            return self.total_cells

        def __getitem__(self, idx):
            fea_1d = self.fea_csr[idx].toarray().flatten().astype(np.float32)
            cell_x = torch.from_numpy(fea_1d).reshape(-1, 1)
            x_fea = torch.from_numpy(fea_1d).unsqueeze(0)
            adj_raw_1d = self.adj_raw_csr[idx].toarray().flatten().astype(np.float32)
            x_adj = torch.from_numpy(adj_raw_1d).unsqueeze(0)
            edge_row = self.adj_edge_csr[idx]
            if edge_row.nnz > 0:
                cols = edge_row.indices
                weights = edge_row.data.astype(np.float32)
                src = self.global_src[cols]
                dst = self.global_dst[cols]
                edge_index = torch.from_numpy(np.stack([src, dst], axis=0)).long()
                edge_attr = torch.from_numpy(weights.reshape(-1, 1))
            else:
                edge_index = torch.empty((2, 0), dtype=torch.long)
                edge_attr = torch.empty((0, 1), dtype=torch.float32)
            y = torch.tensor([self.labels[idx]], dtype=torch.long)
            return Data(x=cell_x, edge_index=edge_index, edge_attr=edge_attr, y=y, x_fea=x_fea, x_adj=x_adj)

    def train_step(svi, train_loader, device):
        epoch_loss = 0.0
        for x_gra in train_loader:
            if 'cuda' in str(device):
                x_gra = x_gra.to(device)
            loss = svi.step(x_gra)
            epoch_loss += loss
        return epoch_loss / max(len(train_loader.dataset), 1)

    def run_train(in_path_gp, in_path_fea, out_path, params, device, pretrain_fea=None, pretrain_adj=None):
        pyro.clear_param_store()
        gc.collect()
        torch.cuda.empty_cache()
        print('=' * 65)
        print('>>> [DOLPHIN Model] High-Efficiency Streaming Engine (RTX 3060 Safe)...')
        print('=' * 65)
        data_dir = os.path.dirname(in_path_fea)
        adj_raw_file = os.path.join(data_dir, params.get('adjacency_raw', 'AdjacencyCompReHvg_fsla.h5ad'))
        adj_edge_file = os.path.join(data_dir, params.get('adjacency_edge', 'AdjacencyCompReHvgEdge_fsla.h5ad'))
        dataset = DolphinLazyDataset(fea_h5ad_path=in_path_fea, adj_raw_h5ad_path=adj_raw_file, adj_edge_h5ad_path=adj_edge_file)
        batch_size = params.get('batch', 4)
        if len(dataset) < 2:
            raise ValueError('Batch-normalized training requires at least two cells')
        batch_size = min(batch_size, len(dataset))
        while len(dataset) % batch_size == 1 and batch_size > 2:
            batch_size -= 1
        if len(dataset) % batch_size == 1:
            raise ValueError('Choose a batch size that avoids a final singleton batch')
        train_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=2, pin_memory=True)
        embedding_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=2)
        in_fea_dim = dataset.num_nodes
        in_adj_dim = dataset.total_adj_cols
        print(f'Dataset Loaded: {dataset.total_cells} cells | Nodes: {in_fea_dim} | Adjacency: {in_adj_dim}')
        print(f'Initializing VAE Model on {device} (Batch size = {batch_size})...')
        vae, svi = define_svi(1, in_fea_dim, in_adj_dim, params, device)
        epochs = params.get('epochs', 100)
        print(f'Starting Training for {epochs} Epochs...')
        epoch_pbar = tqdm(range(epochs), desc='Training Epochs')
        for epoch in epoch_pbar:
            total_loss = train_step(svi, train_loader, device)
            epoch_pbar.set_postfix({'ELBO_Loss': f'{total_loss:.4e}'})
            if (epoch + 1) % 5 == 0:
                gc.collect()
                torch.cuda.empty_cache()
        print('\nExtracting 30-D latent representations (X_z)...')
        vae.eval()
        z_list = []
        with torch.no_grad():
            for x_gra in tqdm(embedding_loader, desc='Inference X_z'):
                if 'cuda' in str(device):
                    x_gra = x_gra.to(device)
                z_mu, _ = vae.getZ(x_gra)
                z_list.append(z_mu.cpu().numpy())
        X_z = np.concatenate(z_list, axis=0)
        adata_z = anndata.AnnData(obs=dataset.adata_fea.obs.copy())
        adata_z.obsm['X_z'] = X_z
        out_file = os.path.join(out_path, 'DOLPHIN_Z.h5ad')
        print(f'Saving latent embedding matrix to: {out_file}')
        adata_z.write_h5ad(out_file, compression='gzip')
        torch.save(vae.state_dict(), os.path.join(out_path, 'model_state.pt'))
        pyro.get_param_store().save(os.path.join(out_path, 'pyro_parameters.pt'))
        print('=' * 65)
        print('>>> DOLPHIN Model Training Finished Successfully!')
        print(f'Output Matrix Shape: {adata_z.shape} | X_z Dimensions: {X_z.shape}')
        print('=' * 65)
    return locals()

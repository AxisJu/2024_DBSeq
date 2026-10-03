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
    import scanpy as sc
    import anndata
    from scipy import sparse
    import warnings
    warnings.filterwarnings('ignore')

    def run_feature_hvg(out_name: str, out_directory: str='./'):
        print('=' * 65)
        print('>>> [Step 8] Normalizing Exon Feature Matrix (ALL Genes Retained)...')
        print('=' * 65)
        final_out_dir = os.path.join(out_directory, 'data')
        feature_anndata = os.path.join(final_out_dir, f'FeatureComp_{out_name}.h5ad')
        print(f'Loading Feature matrix from: {feature_anndata} ...')
        adata = inputs['adata_1'][str(feature_anndata)]
        if not sparse.isspmatrix_csr(adata.X):
            adata.X = sparse.csr_matrix(adata.X)
        print('Applying total count normalization...')
        sc.pp.normalize_total(adata)
        print('Rounding normalized exon counts in sparse format...')
        adata.X.data = np.round(adata.X.data).astype(np.float32)
        out_file = os.path.join(final_out_dir, f'FeatureCompHvg_{out_name}.h5ad')
        print(f'Saving normalized feature matrix to: {out_file} ...')
        adata.write_h5ad(out_file, compression='gzip')
        print(f"Step 8 Finished! Retained {len(set(adata.var['gene_id']))} genes | Shape: {adata.shape}")
    return locals()

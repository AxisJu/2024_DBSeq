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
    import shutil
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

    def _read_cell_batch_worker(batch_info):
        batch_idx, batch_cells, cell_dir = batch_info
        X_list = []
        obs_list = []
        for cell_name in batch_cells:
            cell_file = os.path.join(cell_dir, f'{cell_name}.h5ad')
            if not inputs['data_1']:
                continue
            try:
                ad = inputs['ad_2'][str(cell_file)]
                X_mat = ad.X if sparse.isspmatrix_csr(ad.X) else sparse.csr_matrix(ad.X)
                X_list.append(X_mat)
                obs_list.append(cell_name)
            except Exception as e:
                print(f'Warning: could not read {cell_file}: {e}')
                continue
        if len(X_list) == 0:
            return (batch_idx, None, [])
        batch_X = sparse.vstack(X_list, format='csr')
        return (batch_idx, batch_X, obs_list)

    def run_adjacency_compress_combination(metadata_path: str, out_name: str, out_directory: str='./', adj_run_num: int=50, clean_temp: bool=False, parallel: bool=True, num_processes: int=40):
        print('=' * 65)
        print('>>> [Step 5] Combining Compressed Adjacency Matrices...')
        print('=' * 65)
        final_out_dir = os.path.join(out_directory, 'data')
        cell_dir = os.path.join(final_out_dir, 'temp', 'adj_comp_matrix')
        if not inputs['data_3']:
            raise FileNotFoundError(f'Directory not found: {cell_dir}')
        df_label = inputs['df_label_4'][str(metadata_path)]
        sample_list = df_label['CB'].tolist() if 'CB' in df_label.columns else df_label.iloc[:, 0].tolist()
        total_samples = len(sample_list)
        print(f'Total cells to combine from metadata: {total_samples}')
        print('Extracting feature metadata (var) from first valid cell...')
        var_df = None
        for s in sample_list:
            test_file = os.path.join(cell_dir, f'{s}.h5ad')
            if inputs['data_5'] and os.path.getsize(test_file) > 1024:
                ad_temp = inputs['ad_temp_6'][str(test_file)]
                var_df = ad_temp.var.copy()
                break
        if var_df is None:
            raise RuntimeError('No valid single-cell .h5ad files found in output directory!')
        batches = [(idx, sample_list[i:i + adj_run_num], cell_dir) for idx, i in enumerate(range(0, total_samples, adj_run_num))]
        print(f'Split {total_samples} cells into {len(batches)} batches (Batch size = {adj_run_num})')
        results = []
        workers = num_processes if parallel else 1
        print(f'Streaming batches into memory using {workers} workers...')
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(_read_cell_batch_worker, b) for b in batches]
            for f in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc='Loading Batches'):
                results.append(f.result())
        print('Sorting and executing single-step sparse vstack...')
        results.sort(key=lambda x: x[0])
        valid_X_list = [r[1] for r in results if r[1] is not None]
        all_obs_names = [name for r in results for name in r[2]]
        if len(valid_X_list) == 0:
            raise RuntimeError('No matrix data was loaded. Please check single-cell files.')
        final_X = sparse.vstack(valid_X_list, format='csr')
        final_obs = pd.DataFrame(index=all_obs_names)
        final_adata = anndata.AnnData(X=final_X, obs=final_obs, var=var_df)
        final_output_path = os.path.join(final_out_dir, f'AdjacencyComp_{out_name}.h5ad')
        print(f'Writing final combined AnnData to: {final_output_path}')
        final_adata.write_h5ad(final_output_path, compression='gzip')
        print(f'Successfully generated {final_output_path} | Shape: {final_adata.shape}')
        if clean_temp:
            print('Cleaning up temporary cell directory...')
            if os.path.islink(cell_dir):
                target = inputs['target_7'][str(cell_dir)]
                shutil.rmtree(target, ignore_errors=True)
                os.unlink(cell_dir)
            elif os.path.isdir(cell_dir):
                shutil.rmtree(cell_dir, ignore_errors=True)
        print('Step 5 Finished Successfully!')
    return locals()

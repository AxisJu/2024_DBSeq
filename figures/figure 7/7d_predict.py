def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import os
    from pathlib import Path
    DATA_ROOT = Path(os.environ['DBSEQ_DATA_ROOT'])
    OUTPUT_ROOT = Path(os.environ['DBSEQ_OUTPUT_ROOT'])
    import sys
    import os
    import warnings
    warnings.filterwarnings('ignore')
    import torch
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr, spearmanr
    import scanpy as sc
    matplotlib.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    OUT = str(DATA_ROOT / 'results/260623_squidff')
    DATA = f'{OUT}/data'
    MODEL_DIR = f'{OUT}/models_hvg2000'
    out_dir_fig = str(OUTPUT_ROOT / 'figures/figure 7')
    out_dir_csv = str(OUTPUT_ROOT / 'sourcedata/figure 7')
    os.makedirs(out_dir_fig, exist_ok=True)
    os.makedirs(out_dir_csv, exist_ok=True)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')
    from Squidiff.script_util import model_and_diffusion_defaults, create_model_and_diffusion, args_to_dict
    args = model_and_diffusion_defaults()
    args.update({'gene_size': 2000, 'output_dim': 2000, 'num_layers': 3, 'use_encoder': True, 'diffusion_steps': 1000, 'noise_schedule': 'linear', 'use_drug_structure': False, 'drug_dimension': 1024, 'comb_num': 1})
    model, diffusion = create_model_and_diffusion(**args_to_dict(args, model_and_diffusion_defaults().keys()))
    model.load_state_dict(inputs['data_1'][str(f'{MODEL_DIR}/model.pt')])
    model.to(device)
    model.eval()
    print('Model loaded.')
    train_genes = inputs['train_genes_2'][str(f'{DATA}/atn_dbs_vs_sham_hvg2000.h5ad')].var_names.tolist()
    train_genes = np.array(train_genes)
    adata_full = inputs['adata_full_3'][str(f'{DATA}/atn_training.h5ad')]
    sham_full = inputs['sham_full_4'][str(f'{DATA}/sham_12regions.h5ad')]
    all_genes_f = adata_full.var_names.tolist()
    sham_genes = sham_full.var_names.tolist()

    def align(X, glist):
        g2i = {g: i for i, g in enumerate(glist)}
        out = np.zeros((X.shape[0], 2000), dtype=np.float32)
        for i, g in enumerate(train_genes):
            if g in g2i:
                out[:, i] = X[:, g2i[g]]
        return out

    def ddim_reverse_encode(x0_tensor, z_sem, n_steps=50):
        N = x0_tensor.shape[0]
        z_mod = z_sem.expand(N, -1) if z_sem.shape[0] == 1 else z_sem
        step_size = diffusion.num_timesteps // n_steps
        indices = list(range(0, diffusion.num_timesteps, step_size))
        img = x0_tensor.clone()
        with torch.no_grad():
            for i in indices:
                t = torch.tensor([i] * N, device=device)
                out = diffusion.ddim_reverse_sample(model, img, t, clip_denoised=False, model_kwargs={'z_mod': z_mod, 'x_start': None, 'drug_dose': None, 'group': None})
                img = out['sample']
        return img

    def ddim_decode(z_in, xT=None, n_samples=200):
        if xT is None:
            xT = torch.randn(n_samples, 2000, device=device)
        N = xT.shape[0]
        with torch.no_grad():
            return diffusion.ddim_sample_loop(model, shape=xT.shape, noise=xT, model_kwargs={'z_mod': z_in.expand(N, -1), 'x_start': None, 'drug_dose': None, 'group': None})
    Xa = np.array(adata_full.X).astype(np.float32)
    dbs_mask = adata_full.obs['Group_L2'] == 'DBS_I'
    sham_mask = adata_full.obs['Group_L2'] == 'Sham_I'
    X_dbs_a = torch.tensor(align(Xa[dbs_mask], all_genes_f)).float().to(device)
    X_sham_a = torch.tensor(align(Xa[sham_mask], all_genes_f)).float().to(device)
    real_dbs_mean = X_dbs_a.mean(0).cpu().numpy()
    with torch.no_grad():
        z_dbs = model.encoder(X_dbs_a, label=None, drug_dose=None)
        z_sham = model.encoder(X_sham_a, label=None, drug_dose=None)
    mz_dbs = z_dbs.mean(0)
    mz_sham = z_sham.mean(0)
    dz = mz_dbs - mz_sham
    opt_s = float(os.environ.get('DBSEQ_CALIBRATION_SCALE', '15.0'))
    print(f'Frozen calibration scale: {opt_s}')
    N_rc = min(50, X_sham_a.shape[0])
    X_sham_sub = X_sham_a[:N_rc]
    z_sham_sub = z_sham[:N_rc]
    print(f'Virtual DBS prediction (Sham + dz*{opt_s})...')
    xT_sham = ddim_reverse_encode(X_sham_sub, z_sham_sub, n_steps=20)
    x_virt = ddim_decode(mz_sham.unsqueeze(0) + dz.unsqueeze(0) * opt_s, xT=xT_sham)
    virt_mean = x_virt.mean(0).cpu().numpy()
    virt_sc = (virt_mean - virt_mean.mean()) * real_dbs_mean.std() / (virt_mean.std() + 1e-08) + real_dbs_mean.mean()
    src_df = pd.DataFrame({'gene': train_genes, 'real_dbs_mean': real_dbs_mean, 'virtual_dbs_raw': virt_mean, 'virtual_dbs_scaled': virt_sc})
    src_df.to_csv(f'{out_dir_csv}/7d.csv', index=False)
    import subprocess
    inputs['data_5']
    return locals()

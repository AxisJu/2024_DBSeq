def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']

    def main():
        p = parser(__doc__)
        args = inputs['args_3']
        paths = Paths(args)
        import os
        import scvi
        import torch
        import numpy as np
        import scanpy as sc
        import pandas as pd
        import anndata as ad
        scvi.settings.seed = args.seed
        torch.set_float32_matmul_precision('high')
        adata = inputs['adata_4'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L1.h5ad'))]
        adata.obs['leiden_celltype_level2'] = adata.obs['celltype_level1']
        adata_ExN = inputs['adata_ExN_5'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L1-ExN.h5ad'))]
        exn_map = adata_ExN.obs['leiden_ExN_res_0.3'].astype(str).map(lambda x: f'ExN_{x}')
        colname = 'leiden_celltype_level2'
        existing_categories = adata.obs[colname].cat.categories
        new_categories = pd.Index(exn_map.unique())
        to_add = new_categories.difference(existing_categories)
        adata.obs[colname] = adata.obs[colname].cat.add_categories(to_add)
        adata.obs.loc[adata_ExN.obs_names, colname] = exn_map
        adata_InN = inputs['adata_InN_6'][str(paths.input('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L1-InN.h5ad'))]
        inn_map = adata_InN.obs['leiden_InN_res_0.4'].astype(str).map(lambda x: f'InN_{x}')
        colname = 'leiden_celltype_level2'
        existing_categories = adata.obs[colname].cat.categories
        new_categories = pd.Index(inn_map.unique())
        to_add = new_categories.difference(existing_categories)
        adata.obs[colname] = adata.obs[colname].cat.add_categories(to_add)
        adata.obs.loc[adata_InN.obs_names, colname] = inn_map
        adata.obs['celltype_level2'] = adata.obs['leiden_celltype_level2'].map({'OPC': 'OPC', 'Oligo': 'Oligo', 'Micro': 'Micro', 'Astro': 'Astro', 'VCs': 'VCs', 'EPCs': 'EPCs', 'CHPCs': 'CHPCs', 'ExN_0': 'ExN_TH', 'ExN_1': 'ExN_TH', 'ExN_2': 'ExN_TH', 'ExN_3': 'ExN_CTXpl', 'ExN_4': 'ExN_CTXpl', 'ExN_5': 'ExN_CTXpl', 'ExN_6': 'ExN_TH', 'ExN_7': 'ExN_PAL', 'ExN_8': 'ExN_TH', 'ExN_9': 'ExN_TH', 'ExN_10': 'ExN_CTXpl', 'ExN_11': 'ExN_TH', 'ExN_12': 'ExN_CTXpl', 'ExN_13': 'ExN_TH', 'ExN_14': 'ExN_TH', 'ExN_15': 'ExN_CTXpl', 'ExN_16': 'ExN_PAL', 'ExN_17': 'ExN_CTXsp', 'ExN_18': 'ExN_CTXpl', 'InN_0': 'InN_HY', 'InN_1': 'InN_STR', 'InN_2': 'InN_STR', 'InN_3': 'InN_Nonspecific', 'InN_4': 'InN_STR', 'InN_5': 'InN_TH', 'InN_6': 'InN_PAL', 'InN_7': 'InN_STR', 'InN_8': 'InN_STR', 'InN_9': 'InN_CTXpl', 'InN_10': 'InN_Immature', 'InN_11': 'InN_PAL', 'InN_12': 'ExN_TH', 'InN_13': 'ExN_HY', 'InN_14': 'InN_PAL', 'InN_15': 'InN_STR', 'InN_16': 'InN_HY', 'InN_17': 'InN_MB'})
        adata.obs['celltype_level2'] = adata.obs['celltype_level2'].astype('category')
        adata.obs['celltype_level3'] = adata.obs['leiden_celltype_level2'].map({'OPC': 'OPC', 'Oligo': 'Oligo', 'Micro': 'Micro', 'Astro': 'Astro', 'VCs': 'VCs', 'EPCs': 'EPCs', 'CHPCs': 'CHPCs', 'ExN_0': 'ExN_THns', 'ExN_1': 'ExN_ATN', 'ExN_2': 'ExN_CM', 'ExN_3': 'ExN_RSP_L45IT', 'ExN_4': 'ExN_RSP_L6CT', 'ExN_5': 'ExN_DG', 'ExN_6': 'ExN_MH', 'ExN_7': 'ExN_TRS', 'ExN_8': 'ExN_RE', 'ExN_9': 'ExN_ATN', 'ExN_10': 'ExN_CA3', 'ExN_11': 'ExN_LH', 'ExN_12': 'ExN_CA1', 'ExN_13': 'ExN_PF', 'ExN_14': 'ExN_ATN', 'ExN_15': 'ExN_RSP_L23IT', 'ExN_16': 'ExN_BST', 'ExN_17': 'ExN_CLA', 'ExN_18': 'ExN_HPF_CajalRetzius', 'InN_0': 'InN_HYa', 'InN_1': 'InN_CP_D1', 'InN_2': 'InN_CP_D2', 'InN_3': 'InN_Nonspecific', 'InN_4': 'InN_LSc', 'InN_5': 'InN_RT', 'InN_6': 'InN_GPe_MGE', 'InN_7': 'InN_STRns_LGE', 'InN_8': 'InN_STRns_MGE', 'InN_9': 'InN_RSP_MGE', 'InN_10': 'InN_Immature', 'InN_11': 'InN_GPi_MGE', 'InN_12': 'ExN_ATN', 'InN_13': 'ExN_HYa', 'InN_14': 'InN_GPe_LGE', 'InN_15': 'InN_OT', 'InN_16': 'InN_HYa', 'InN_17': 'InN_SCsg'})
        adata.obs['celltype_level3'] = adata.obs['celltype_level3'].astype('category')
        adata.obs['brainregion_full'] = adata.obs['leiden_celltype_level2'].map({'OPC': '/', 'Oligo': '/', 'Micro': '/', 'Astro': '/', 'VCs': '/', 'EPCs': '/', 'CHPCs': '/', 'ExN_0': 'TH', 'ExN_1': 'TH_DORpm_ATN_AV', 'ExN_2': 'TH_DORpm_ILM_CM', 'ExN_3': 'CTXpl_IsoCortex_RSP_RSPd', 'ExN_4': 'CTXpl_IsoCortex_RSP_RSPd', 'ExN_5': 'CTXpl_HPF_HIP_DG', 'ExN_6': 'TH_DORpm_EPI_MH', 'ExN_7': 'PAL_PALc_TRS', 'ExN_8': 'TH_DORpm_RE', 'ExN_9': 'TH_DORpm_ATN_AV', 'ExN_10': 'CTXpl_HPF_HIP_CA_CA3', 'ExN_11': 'TH_DORpm_EPI_LH', 'ExN_12': 'CTXpl_HPF_HIP_CA_CA1', 'ExN_13': 'TH_DORpm_ILM_PF', 'ExN_14': 'TH_DORpm_ATN_AD', 'ExN_15': 'CTXpl_IsoCortex_RSP_RSPv', 'ExN_16': 'PAL_PALc_BST', 'ExN_17': 'CTXsp_CLA', 'ExN_18': 'CTXpl_HPF', 'InN_0': 'HY', 'InN_1': 'STR_STRd_CP', 'InN_2': 'STR_STRd_CP', 'InN_3': '/', 'InN_4': 'STR_LSX_LS_LSc', 'InN_5': 'TH_DORpm_RT', 'InN_6': 'PAL_PALd_GPe', 'InN_7': 'STR', 'InN_8': 'STR', 'InN_9': 'CTXpl_IsoCortex_RSP_RSPd', 'InN_10': '/', 'InN_11': 'PAL_PALd_GPi', 'InN_12': 'TH_DORpm_ATN_AV', 'InN_13': 'HY', 'InN_14': 'PAL_PALd_Gpe', 'InN_15': 'STR_STRv_OT', 'InN_16': 'HY', 'InN_17': 'MB_MBsen_SCs_SCsg'})
        adata.obs['brainregion_full'] = adata.obs['brainregion_full'].astype('category')
        adata.obs['brainregion_level3'] = adata.obs['leiden_celltype_level2'].map({'OPC': '/', 'Oligo': '/', 'Micro': '/', 'Astro': '/', 'VCs': '/', 'EPCs': '/', 'CHPCs': '/', 'ExN_0': 'TH', 'ExN_1': 'TH', 'ExN_2': 'TH', 'ExN_3': 'CTXpl', 'ExN_4': 'CTXpl', 'ExN_5': 'CTXpl', 'ExN_6': 'TH', 'ExN_7': 'PAL', 'ExN_8': 'TH', 'ExN_9': 'TH', 'ExN_10': 'CTXpl', 'ExN_11': 'TH', 'ExN_12': 'CTXpl', 'ExN_13': 'TH', 'ExN_14': 'TH', 'ExN_15': 'CTXpl', 'ExN_16': 'PAL', 'ExN_17': 'CTXsp', 'ExN_18': 'CTXpl', 'InN_0': 'HY', 'InN_1': 'STR', 'InN_2': 'STR', 'InN_3': '/', 'InN_4': 'STR', 'InN_5': 'TH', 'InN_6': 'PAL', 'InN_7': 'STR', 'InN_8': 'STR', 'InN_9': 'CTXpl', 'InN_10': '/', 'InN_11': 'PAL', 'InN_12': 'TH', 'InN_13': 'HY', 'InN_14': 'PAL', 'InN_15': 'STR', 'InN_16': 'HY', 'InN_17': 'MB'})
        adata.obs['brainregion_level3'] = adata.obs['brainregion_level3'].astype('category')
        adata.obs['brainregion_level4'] = adata.obs['leiden_celltype_level2'].map({'OPC': '/', 'Oligo': '/', 'Micro': '/', 'Astro': '/', 'VCs': '/', 'EPCs': '/', 'CHPCs': '/', 'ExN_0': '/', 'ExN_1': 'DORpm', 'ExN_2': 'DORpm', 'ExN_3': 'IsoCortex', 'ExN_4': 'IsoCortex', 'ExN_5': 'HPF', 'ExN_6': 'DORpm', 'ExN_7': 'PALc', 'ExN_8': 'DORpm', 'ExN_9': 'DORpm', 'ExN_10': 'HPF', 'ExN_11': 'DORpm', 'ExN_12': 'HPF', 'ExN_13': 'DORpm', 'ExN_14': 'DORpm', 'ExN_15': 'IsoCortex', 'ExN_16': 'PALc', 'ExN_17': 'CLA', 'ExN_18': 'HPF', 'InN_0': '/', 'InN_1': 'STRd', 'InN_2': 'STRd', 'InN_3': '/', 'InN_4': 'LSX', 'InN_5': 'DORpm', 'InN_6': 'PALd', 'InN_7': '/', 'InN_8': '/', 'InN_9': 'IsoCortex', 'InN_10': '/', 'InN_11': 'PALd', 'InN_12': 'DORpm', 'InN_13': '/', 'InN_14': 'PALd', 'InN_15': 'STRv', 'InN_16': '/', 'InN_17': 'MBsen'})
        adata.obs['brainregion_level4'] = adata.obs['brainregion_level4'].astype('category')
        adata.obs['brainregion_projection'] = adata.obs['leiden_celltype_level2'].map({'OPC': '/', 'Oligo': '/', 'Micro': '/', 'Astro': '/', 'VCs': '/', 'EPCs': '/', 'CHPCs': '/', 'ExN_0': '/', 'ExN_1': 'AV', 'ExN_2': 'CM', 'ExN_3': 'RSPd', 'ExN_4': 'RSPd', 'ExN_5': 'DG', 'ExN_6': 'MH', 'ExN_7': 'TRS', 'ExN_8': 'RE', 'ExN_9': 'AV', 'ExN_10': 'CA3', 'ExN_11': 'LH', 'ExN_12': 'CA1', 'ExN_13': 'PF', 'ExN_14': 'AD', 'ExN_15': 'RSPv', 'ExN_16': 'BST', 'ExN_17': 'CLA', 'ExN_18': '/', 'InN_0': '/', 'InN_1': 'CP', 'InN_2': 'CP', 'InN_3': '/', 'InN_4': 'LSc', 'InN_5': 'RT', 'InN_6': 'GPe', 'InN_7': '/', 'InN_8': '/', 'InN_9': 'RSPd', 'InN_10': '/', 'InN_11': 'GPi', 'InN_12': 'AV', 'InN_13': '/', 'InN_14': 'GPe', 'InN_15': 'OT', 'InN_16': '/', 'InN_17': 'SCs'})
        adata.obs['brainregion_projection'] = adata.obs['brainregion_projection'].astype('category')
        adata.write(paths.output('data/snRNAseq_mouse/processed/matrix/merged_preprocessed_scvi_umap_leiden_annotation_L23.h5ad'), compression='gzip')
        adata.obs.to_csv(paths.output('data/snRNAseq_mouse/processed/metadata/metadata_obs_107.csv'), index=True)
    main()
    return locals()

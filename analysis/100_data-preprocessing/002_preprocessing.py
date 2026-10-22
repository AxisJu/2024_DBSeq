def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']
    SAMPLES = [{'Sample': 'DBS_I_1', 'Group_L1': 'DBS', 'Group_L2': 'DBS_I', 'Sample_Time': '20240805', 'Sequencing_Time': '20240812', 'folder': 'DBS-1', 'max_counts': 10000, 'max_genes': 5000}, {'Sample': 'DBS_I_2', 'Group_L1': 'DBS', 'Group_L2': 'DBS_I', 'Sample_Time': '20241109', 'Sequencing_Time': '20241120', 'folder': 'DBS-2', 'max_counts': 50000, 'max_genes': 8000}, {'Sample': 'DBS_C_1', 'Group_L1': 'Control', 'Group_L2': 'DBS_C', 'Sample_Time': '20240805', 'Sequencing_Time': '20240812', 'folder': 'DBSDC-1', 'max_counts': 35000, 'max_genes': 8000}, {'Sample': 'PTZ_1', 'Group_L1': 'Control', 'Group_L2': 'PTZ', 'Sample_Time': '20240423', 'Sequencing_Time': '20240902', 'folder': 'PTZ-old', 'max_counts': 20000, 'max_genes': 6000}, {'Sample': 'Saline_I_1', 'Group_L1': 'Control', 'Group_L2': 'Saline_I', 'Sample_Time': '20241109', 'Sequencing_Time': '20241120', 'folder': 'Saline-1', 'max_counts': 60000, 'max_genes': 9000}, {'Sample': 'Saline_I_2', 'Group_L1': 'Control', 'Group_L2': 'Saline_I', 'Sample_Time': '20241109', 'Sequencing_Time': '20241120', 'folder': 'Saline-2', 'max_counts': 70000, 'max_genes': 10000}, {'Sample': 'Sham_I_1', 'Group_L1': 'Control', 'Group_L2': 'Sham_I', 'Sample_Time': '20241219', 'Sequencing_Time': '20241224', 'folder': 'Sham-1', 'max_counts': 40000, 'max_genes': 8000}, {'Sample': 'Sham_I_2', 'Group_L1': 'Control', 'Group_L2': 'Sham_I', 'Sample_Time': '20241219', 'Sequencing_Time': '20241224', 'folder': 'Sham-2', 'max_counts': 60000, 'max_genes': 8000}, {'Sample': 'Sham_I_3', 'Group_L1': 'Control', 'Group_L2': 'Sham_I', 'Sample_Time': '20241219', 'Sequencing_Time': '20241224', 'folder': 'Sham-3', 'max_counts': 60000, 'max_genes': 8000}]

    def main():
        p = parser(__doc__)
        args = inputs['args_3']
        paths = Paths(args)
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import scrublet as scr
        from scipy import sparse
        np.random.seed(args.seed)
        base = 'data/snRNAseq_mouse/processed'
        selected = SAMPLES if args.samples is None else [s for s in SAMPLES if s['Sample'] in args.samples]
        if not selected or (args.samples and len(selected) != len(set(args.samples))):
            raise ValueError('Unknown sample ID')
        objects, qc = ({}, [])
        for spec in selected:
            raw = paths.input(f"data/snRNAseq_mouse/raw/{spec['folder']}/cellranger/filtered_feature_bc_matrix.h5")
            obj = inputs['obj_4'][str(raw)]
            obj.var_names_make_unique()
            for key in ('Sample', 'Group_L1', 'Group_L2', 'Sample_Time', 'Sequencing_Time'):
                obj.obs[key] = spec[key]
            obj.var['mt'] = obj.var_names.str.startswith('mt-')
            obj.var['ribo'] = obj.var_names.str.startswith(('Rps', 'Rpl'))
            obj.var['hb'] = obj.var_names.str.contains('^Hb[^(p)]')
            sc.pp.calculate_qc_metrics(obj, qc_vars=['mt', 'ribo', 'hb'], inplace=True, log1p=True)
            scrub = scr.Scrublet(obj.X, random_state=args.seed)
            scores, doublets = scrub.scrub_doublets()
            obj.obs['doublet_scores'] = scores
            obj.obs['predicted_doublets'] = doublets
            initial = obj.n_obs
            sc.pp.filter_cells(obj, min_genes=200)
            sc.pp.filter_genes(obj, min_cells=5)
            keep = (obj.obs.total_counts < spec['max_counts']) & (obj.obs.n_genes_by_counts < spec['max_genes']) & (obj.obs.pct_counts_mt < 10) & (obj.obs.pct_counts_ribo < 10) & ~obj.obs.predicted_doublets
            obj = obj[keep].copy()
            obj.layers['counts'] = obj.X.copy()
            sc.pp.normalize_total(obj, target_sum=10000.0)
            sc.pp.log1p(obj)
            obj.layers['data'] = obj.X.copy()
            objects[spec['Sample']] = obj
            qc.append(dict(Sample=spec['Sample'], cells_before=initial, cells_after=obj.n_obs))
        merged = sc.concat(objects, join='inner', label='Sample')
        raw_names = merged.obs_names.astype(str).copy()
        merged.obs_names_make_unique()
        pd.DataFrame(dict(raw_cellname=raw_names, unique_cellname=merged.obs_names)).to_csv(paths.output(base + '/metadata/metadata_cellnamemapping.csv'), index=False)
        merged.raw = merged.copy()
        merged.write_h5ad(paths.output(base + '/matrix/merged_preprocessed.h5ad'), compression='gzip')
        for layer, filename in [('counts', 'counts'), ('data', 'logdata')]:
            sparse.save_npz(paths.output(base + f'/matrix/matrix_{filename}.npz'), sparse.csc_matrix(merged.layers[layer].T))
        merged.obs.to_csv(paths.output(base + '/metadata/metadata_obs_101.csv'))
        merged.var.to_csv(paths.output(base + '/metadata/metadata_var.csv'))
        pd.DataFrame(qc).to_csv(paths.output('results/Table/preprocessing_qc.csv'), index=False)
    main()
    return locals()

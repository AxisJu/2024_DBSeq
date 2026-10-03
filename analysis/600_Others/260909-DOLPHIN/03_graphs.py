def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']
    load_adaptation = inputs['load_adaptation_3']
    patch_original = inputs['patch_original_4']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_5']
        data, out = roots('600_Others/260909-DOLPHIN')
        import pandas as pd
        import sys
        if sys.platform == 'win32' and a.stage in [1, 3, 4, 5]:
            raise RuntimeError('Parallel DOLPHIN graph adaptations require Linux fork workers')
        meta = input_path(data, a.metadata)
        frame = inputs['frame_6'][str(meta)]
        if not {'CB', 'donor_id', 'Group_L2', 'final_region'} <= set(frame.columns) or frame.CB.duplicated().any():
            raise ValueError('Unique cell metadata with biological donor IDs is required')
        ref = input_path(data, a.reference)
        common = dict(out_directory=str(out), out_name='fsla')
        if a.stage == 1:
            from DOLPHIN.graph_generation.preprocess_raw_reads import run_parallel_gene_processing
            run_parallel_gene_processing(metadata_path=str(meta), gtf_path=str(ref / 'dolphin.exon.pkl'), adj_index_path=str(ref / 'dolphin_adj_index.csv'), main_folder=str(out), n_processes=a.workers)
        elif a.stage == 2:
            inputs['data_7']
            from DOLPHIN.graph_generation.process_feature_matrix import run_feature_combination
            run_feature_combination(metadata_path=str(meta), graph_directory=str(out / '06_graph_mtx'), gene_annotation=str(ref / 'dolphin_gene_meta.csv'), gtf_pkl_path=str(ref / 'dolphin.exon.pkl'), clean_temp=False, **common)
        elif a.stage == 3:
            inputs['data_8']
            module = inputs['module_9']
            module.run_adjacency_combination(metadata_path=str(meta), graph_directory=str(out / '06_graph_mtx'), adj_meta_file=str(ref / 'dolphin_adj_metadata_table.csv'), clean_temp=False, adj_run_num=50, parallel=False, **common)
        elif a.stage == 4:
            inputs['data_10'][str('process_adjacency_matrix_compress')].run_adjacency_compression(metadata_path=str(meta), num_processes=a.workers, **common)
        elif a.stage == 5:
            inputs['data_11'][str('process_adjacency_matrix_compress_combine')].run_adjacency_compress_combination(metadata_path=str(meta), num_processes=a.workers, clean_temp=False, **common)
        elif a.stage == 6:
            inputs['data_12'][str('process_adjacency_matrix_final')].run_adjacency_matrix_final(**common)
        elif a.stage == 7:
            inputs['data_13'][str('process_feature_hvg')].run_feature_hvg(**common)
        elif a.stage == 8:
            inputs['data_14'][str('process_adjacency_hvg')].run_adjacency_hvg(**common)
        elif a.stage == 9:
            inputs['data_15'][str('process_graph_final')].run_model_input(metadata_path=str(meta), celltypename='final_region', **common)
    main()
    return locals()

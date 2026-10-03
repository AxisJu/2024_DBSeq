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
        import subprocess
        import pandas as pd
        base = 'data/snRNAseq_mouse/processed/intermediate/scDRS'
        magma = 'results/MAGMA'
        counts = base + '/001_input/running_all_250704_onlycounts_obs.h5ad'
        covariates = base + '/001_input/running_all_250704_onlycounts.cov.tsv'
        gs = base + f'/001_input/gs_file/{args.trait}.gs'
        scores = magma + f'/step2/{args.trait}.tsv'
        command = None
        if args.stage == 'annotate':
            command = [args.magma, '--annotate', 'window=10,10', '--snp-loc', str(paths.input(args.bfile + '.bim')), '--gene-loc', str(paths.input(args.gene_locations)), '--out', str(paths.output(magma + '/step1'))]
        elif args.stage == 'magma':
            if not args.input:
                p.error('magma requires --input GWAS summary statistics')
            command = [args.magma, '--bfile', str(paths.input(args.bfile)), '--pval', str(paths.input(args.input)), 'use=MarkerName,P-value', 'ncol=Effective_N', '--gene-annot', str(paths.input(magma + '/step1.genes.annot')), '--out', str(paths.output(magma + f'/step2/{args.trait}'))]
        elif args.stage == 'gene-scores':
            locations = inputs['locations_4'][str(paths.input(args.gene_locations))]
            stats = inputs['stats_5'][str(paths.input(args.input or magma + f'/step2/{args.trait}.genes.out'))]
            result = stats.merge(locations[['ENTREZ', 'SYMBOL']], left_on='GENE', right_on='ENTREZ', how='left', validate='many_to_one')
            result = result[['SYMBOL', 'ZSTAT']].dropna().rename(columns={'SYMBOL': 'GENE', 'ZSTAT': args.trait + '_Z'})
            if result.GENE.duplicated().any():
                raise ValueError('Gene symbols map to multiple MAGMA entries; resolve the reference mapping')
            result.to_csv(paths.output(scores), sep='\t', index=False)
        elif args.stage == 'munge':
            command = [args.scdrs, 'munge-gs', '--out-file', str(paths.output(gs)), '--zscore-file', str(paths.input(scores)), '--weight', 'zscore', '--fdr', '0.05', '--n-max', '1000']
        elif args.stage == 'prepare':
            import anndata
            obj = inputs['obj_6'][str(paths.input(args.input or 'data/snRNAseq_mouse/processed/matrix/running_all_250704.h5ad'))]
            raw = anndata.AnnData(obj.layers['counts'].copy(), obs=obj.obs.copy(), var=obj.var.copy())
            raw.write_h5ad(paths.output(counts))
            cov = obj.obs[['n_genes_by_counts', 'total_counts']].rename_axis('index')
            cov.to_csv(paths.output(covariates), sep='\t')
        elif args.stage == 'score':
            command = [args.scdrs, 'compute-score', '--h5ad-file', str(paths.input(counts)), '--h5ad-species', 'mmusculus', '--gs-file', str(paths.input(gs)), '--gs-species', 'hsapiens', '--cov-file', str(paths.input(covariates)), '--flag-filter-data', 'True', '--flag-raw-count', 'True', '--flag-return-ctrl-raw-score', 'False', '--flag-return-ctrl-norm-score', 'True', '--out-folder', str(paths.output(base + '/002_output'))]
        elif args.stage == 'downstream':
            command = [args.scdrs, 'perform-downstream', '--h5ad-file', str(paths.input(counts)), '--score-file', str(paths.input(base + f'/002_output/{args.trait}_Z.full_score.gz')), '--out-folder', str(paths.output(base + '/002_output')), '--group-analysis', args.group_column, '--flag-filter-data', 'True', '--flag-raw-count', 'True']
        if command:
            print(subprocess.list2cmdline(command))
            if not args.dry_run:
                inputs['data_7']
    main()
    return locals()

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']
    LIBRARIES = {'DBS-1': 3, 'DBS-2': 2, 'DBSDC-1': 3, 'PTZ-old': 1, 'Saline-1': 3, 'Saline-2': 3, 'Sham-1': 3, 'Sham-2': 3, 'Sham-3': 3}

    def main():
        p = parser(__doc__)
        args = inputs['args_3']
        paths = Paths(args)
        import subprocess
        base = 'data/snRNAseq_mouse'
        intermediate = base + '/processed/intermediate'
        bams = [str(paths.input(f'{base}/raw/{lib}/bam/possorted_genome_bam.bam')) for lib in args.libraries]
        commands = []
        if args.stage == 'pooled-snps':
            dest = paths.output(intermediate + '/cellSNP/Model_2b_250223')
            commands.append([args.cellsnp, '-s', ','.join(bams), '-O', str(dest), '-p', str(args.threads), '--minMAF', '0.01', '--minCOUNT', '100', '--cellTAG', 'None', '--UMItag', 'UB', '--chrom', ','.join(map(str, range(1, 20))), '--gzip'])
        for lib, bam in zip(args.libraries, bams):
            cell_dir = intermediate + f'/cellSNP/Model_1a_250228/{lib}'
            if args.stage == 'cell-snps':
                commands.append([args.cellsnp, '-s', bam, '-b', str(paths.input(f'{base}/raw/{lib}/cellranger/filtered_feature_bc_matrix/barcodes.tsv.gz')), '-O', str(paths.output(cell_dir)), '-R', str(paths.input(intermediate + '/cellSNP/Model_2b_250223/cellSNP.base.vcf.gz')), '-p', str(args.threads), '--minMAF', '0.01', '--minCOUNT', '20', '--gzip'])
            if args.stage == 'vireo':
                commands.append([args.vireo, '-c', str(paths.input(cell_dir)), '-N', str(LIBRARIES[lib]), '-o', str(paths.output(intermediate + f'/vireoSNP/{lib}')), '--nInit=50', '--extraDonor=0', f'--randSeed={args.seed}', '--callAmbientRNAs', f'--nproc={args.threads}'])
            if args.stage == 'donor-number' and (not args.dry_run):
                import pandas as pd
                import vireoSNP
                from scipy.io import mmread
                ad = inputs['ad_4'][str(paths.input(cell_dir + '/cellSNP.tag.AD.mtx'))].tocsc()
                dp = inputs['dp_5'][str(paths.input(cell_dir + '/cellSNP.tag.DP.mtx'))].tocsc()
                rows = []
                for n in range(2, 7):
                    result = vireoSNP.vireo_wrap(ad, dp, n_donor=n, learn_GT=True, n_extra_donor=0, ASE_mode=False, fix_beta_sum=False, n_init=50, check_doublet=True, random_seed=args.seed)
                    rows.extend((dict(n_donors=n, initialization=i, elbo=value) for i, value in enumerate(result['LB_list'])))
                pd.DataFrame(rows).to_csv(paths.output(f'results/Table/vireo/{lib}_donor_number.csv'), index=False)
        for command in commands:
            print(subprocess.list2cmdline(command))
            if not args.dry_run:
                inputs['data_6']
    main()
    return locals()

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
        import pandas as pd
        import vireoSNP
        result = vireoSNP.vcf.match_VCF_samples(str(paths.input(args.first)), str(paths.input(args.second)), GT_tag1=args.tag, GT_tag2=args.tag)
        values = pd.DataFrame(result['matched_GPb_diff'], index=result['matched_donors1'], columns=result['matched_donors2'])
        values.to_csv(paths.output(f'results/Table/vireo/{args.name}_distance.csv'))
        pd.DataFrame({'shared_variants': [result['matched_n_var']]}).to_csv(paths.output(f'results/Table/vireo/{args.name}_variants.csv'), index=False)
    main()
    return locals()

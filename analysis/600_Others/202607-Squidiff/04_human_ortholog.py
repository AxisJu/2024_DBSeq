def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/202607-Squidiff')
        import pandas as pd
        from scipy.stats import spearmanr
        mapping = inputs['mapping_4'][str(input_path(data, a.orthologs))].drop_duplicates(['mouse_gene', 'human_gene'])
        ambiguous = mapping.mouse_gene.duplicated(False) | mapping.human_gene.duplicated(False)
        mapping.loc[ambiguous].to_csv(out / 'ambiguous_orthologs.csv', index=False)
        mapping = mapping.loc[~ambiguous]
        mouse = inputs['mouse_5'][str(input_path(data, a.mouse))]
        human = inputs['human_6'][str(input_path(data, a.human))]
        if mouse.duplicated(['region', 'gene']).any() or human.duplicated(['cohort', 'human_gene']).any():
            raise ValueError('Duplicate gene effects within a region or cohort')
        aligned = mouse.merge(mapping, left_on='gene', right_on='mouse_gene', validate='many_to_one').merge(human, on='human_gene')
        aligned.to_csv(out / 'ortholog_aligned_effects.csv', index=False)
        rows = []
        for (region, cohort), frame in aligned.groupby(['region', 'cohort']):
            frame = frame.dropna(subset=['encoder_gradient', 'log2FC'])
            if len(frame) >= 3:
                rho, pv = spearmanr(frame.encoder_gradient, frame.log2FC)
                rows.append(dict(region=region, cohort=cohort, n_genes=len(frame), rho=rho, descriptive_gene_p=pv))
        pd.DataFrame(rows).to_csv(out / 'cross_species_descriptive_correlations.csv', index=False)
    main()
    return locals()

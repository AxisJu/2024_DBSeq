def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('500_GRN')
        import pandas as pd
        genes = inputs['genes_4'][str(input_path(data, a.genes))].gene.drop_duplicates()
        if a.edges:
            if not a.mapping:
                p.error('--mapping is required with --edges')
            import networkx as nx
            mapping = inputs['mapping_5'][str(input_path(data, a.mapping))].drop_duplicates(['gene', 'protein_id'])
            links = inputs['links_6'][str(input_path(data, a.edges))]
            links = links[links.combined_score >= a.score]
            graph = nx.from_pandas_edgelist(links, 'protein1', 'protein2')
            rows = []
            for gene in genes:
                proteins = set(mapping.loc[mapping.gene == gene, 'protein_id'])
                neighbors = set().union(*(set(graph.neighbors(k)) for k in proteins if k in graph)) if proteins else set()
                rows.append(dict(gene=gene, mapped=bool(proteins), n_proteins=len(proteins), degree_union=len(neighbors - proteins) if proteins else None, protein_ids=';'.join(sorted(proteins))))
            pd.DataFrame(rows).to_csv(out / 'STRING_degrees.csv', index=False)
        if a.gmt:
            rows = []
            universe = set(genes)
            handle = inputs['handle_7'][str(input_path(data, a.gmt))]
            for line in handle:
                term, description, *members = line.rstrip().split('\t')
                rows.extend(({'gene': gene, 'term': term} for gene in sorted(set(members) & universe)))
            memberships = pd.DataFrame(rows, columns=['gene', 'term']).drop_duplicates()
            memberships.to_csv(out / 'GO_memberships.csv', index=False)
            counts = memberships.groupby('gene').size().reindex(genes, fill_value=0)
            counts.rename('n_GO_terms').rename_axis('gene').to_csv(out / 'GO_membership_counts.csv')
    main()
    return locals()

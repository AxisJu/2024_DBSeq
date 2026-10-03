def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260909-DOLPHIN')
        out = out / 'reference'
        out.mkdir(exist_ok=True)
        new = out / 'genes_dolphin.gtf'
        src = inputs['src_4'][str(input_path(data, a.gtf))]
        dst = inputs['dst_5'][str(new)]
        for line in src:
            if line.startswith('#'):
                dst.write(line)
                continue
            parts = line.rstrip().split('\t')
            if len(parts) != 9:
                raise ValueError('Invalid GTF row')
            parts[8] = parts[8].replace('gene_type "', 'gene_biotype "')
            if 'gene_source "' not in parts[8]:
                parts[8] = parts[8].rstrip(';') + '; gene_source "10X";'
            dst.write('\t'.join(parts) + '\n')
        index = inputs['adjacency_index']
        adj, genes = inputs['adjacency_metadata']
        index.to_csv(out / 'dolphin_adj_index.csv', index=False)
        adj.to_csv(out / 'dolphin_adj_metadata_table.csv', index=False)
        genes.to_csv(out / 'dolphin_gene_meta.csv', index=False)
    main()
    return locals()

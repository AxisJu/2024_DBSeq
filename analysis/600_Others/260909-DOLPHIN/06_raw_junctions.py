def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import csv
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260909-DOLPHIN')
        handle = inputs['handle_4'][str(input_path(data, a.manifest))]
        manifest = list(csv.DictReader(handle))
        if len({r['sample_id'] for r in manifest}) != len(manifest):
            raise ValueError('Duplicate cell entries')
        dst = inputs['dst_5'][str(out / 'raw_junction_reads.csv')]
        writer = csv.DictWriter(dst, fieldnames=['sample_id', 'junction_id', 'reads', 'read_origin'])
        writer.writeheader()
        for row in manifest:
            path = input_path(data, row['sj_path'])
            if '.aggr' in row['sample_id'] or '.aggr' in str(path):
                raise ValueError('Aggregated reads cannot support donor-level inference')
            source = inputs['source_6'][str(path)]
            for line in source:
                fields = line.rstrip().split('\t')
                chrom, start, end, strand = fields[:4]
                if strand not in ['1', '2']:
                    continue
                sign = '+' if strand == '1' else '-'
                writer.writerow(dict(sample_id=row['sample_id'], junction_id=f'junction:{chrom}:{start}-{end}:{sign}', reads=int(fields[6]), read_origin='original_STAR_unique'))
    main()
    return locals()

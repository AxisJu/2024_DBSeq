def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import math
    from collections import Counter
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def read_junctions(read):
        position = read.reference_start
        result = []
        for operation, length in read.cigartuples or []:
            if operation == 3:
                result.append((read.reference_name, position + 1, position + length))
            if operation in [0, 2, 3, 7, 8]:
                position += length
        return result

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260909-DOLPHIN')
        import pandas as pd
        import pysam
        if not 0 < a.junction_frequency <= 1:
            raise ValueError('Frequency must lie in (0,1]')
        meta = inputs['meta_4'][str(input_path(data, a.metadata))].set_index('CB')
        manifest = inputs['manifest_5'][str(input_path(data, a.manifest))].set_index('sample_id')
        neighbors = inputs['neighbors_6'][str(input_path(data, a.neighbors))]
        if not meta.index.is_unique or not manifest.index.is_unique:
            raise ValueError('Duplicate cell IDs')
        destination = out / 'exploratory_smoothed_bams'
        destination.mkdir(exist_ok=True)
        audit = []
        keys = ['donor_id', 'Group_L2', 'final_region']
        for target, frame in neighbors.groupby('CB'):
            related = sorted(set(frame.neighbor) - {target})
            for cell in related:
                if not meta.loc[target, keys].equals(meta.loc[cell, keys]):
                    raise ValueError('Neighbor crosses a donor, condition, or region boundary')
            support = Counter()
            for cell in related:
                handle = inputs['handle_7'][str(input_path(data, manifest.loc[cell, 'sj_path']))]
                junctions = {(r[0], int(r[1]), int(r[2])) for r in (line.rstrip().split('\t') for line in handle)}
                support.update(junctions)
            threshold = max(2, math.ceil(len(related) * a.junction_frequency))
            accepted = {j for j, n in support.items() if n >= threshold}
            unsorted = destination / (str(target) + '.aggr.unsorted.bam')
            target_bam = input_path(data, manifest.loc[target, 'bam'])
            template = inputs['template_8'][str(str(target_bam))]
            dst = inputs['dst_9'][str(str(unsorted))]
            for read in inputs['data_10']:
                dst.write(read)
            for cell in related:
                count = 0
                src = inputs['src_11'][str(str(input_path(data, manifest.loc[cell, 'bam'])))]
                if src.references != template.references:
                    raise ValueError('BAM reference dictionaries differ')
                for read in inputs['data_12']:
                    if set(read_junctions(read)) & accepted:
                        dst.write(read)
                        count += 1
                audit.append(dict(target=target, source_cell=cell, borrowed_reads=count, use='exploratory_only'))
            final = destination / (str(target) + '.aggr.bam')
            inputs['data_13']
            inputs['data_14']
        pd.DataFrame(audit).to_csv(out / 'exploratory_read_borrowing.csv', index=False)
    main()
    return locals()

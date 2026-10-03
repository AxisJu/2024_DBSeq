def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import csv
    import gzip
    import subprocess
    from pathlib import Path
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def read_bed(path, convention):
        op = gzip.open if str(path).endswith('.gz') else open
        with op(path, 'rt') as handle:
            for index, line in enumerate(handle):
                if line.startswith(('#', 'track', 'browser')) or not line.strip():
                    continue
                fields = line.split()
                start, end = (int(fields[1]), int(fields[2]))
                if convention == 'closed1':
                    start -= 1
                if start < 0 or end <= start:
                    raise ValueError('Invalid peak interval')
                yield (fields[0], start, end, index)

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260921-CUT&TAG')
        import pandas as pd
        manifest = inputs['manifest_4'][str(input_path(data, a.manifest))].fillna('')
        if manifest.sample_id.duplicated().any():
            raise ValueError('Duplicate sample IDs')
        rows = []
        for sample in manifest.itertuples():
            if sample.coordinate_system not in ['bed0', 'closed1']:
                raise ValueError('Explicit bed0 or closed1 convention is required')
            intervals = list(inputs['intervals_5'][str(input_path(data, sample.path))])
            if sample.genome != 'mm10':
                if not sample.chain:
                    raise ValueError('A genome-specific liftOver chain is required for ' + sample.genome)
                name = str(sample.sample_id)
                if Path(name).name != name:
                    raise ValueError('Invalid sample ID')
                source, lifted, unmapped = [out / (name + suffix) for suffix in ['.input.bed', '.mm10.bed', '.unmapped.bed']]
                handle = inputs['handle_6'][str(source)]
                for chrom, start, end, identifier in intervals:
                    handle.write(f'{chrom}\t{start}\t{end}\t{identifier}\n')
                inputs['data_7']
                mapped = [line.split() for line in inputs['mapped_8'][str(lifted)].splitlines() if line.strip()]
                from collections import Counter
                multiplicity = Counter((row[3] for row in mapped))
                intervals = [(r[0], int(r[1]), int(r[2]), int(r[3])) for r in mapped if multiplicity[r[3]] == 1]
                pd.DataFrame([{'source_peak': k, 'n_mappings': v} for k, v in multiplicity.items()]).to_csv(out / (name + '.mapping_status.csv'), index=False)
            rows.extend((dict(sample_id=sample.sample_id, assay=sample.assay, genome='mm10', chrom=c, start=s + 1, end=e, source_peak=i, source_genome=sample.genome) for c, s, e, i in intervals))
        original = pd.DataFrame(rows)
        original.to_csv(out / 'mm10_peaks_by_sample.csv', index=False)
        union = []
        for (assay, chrom), frame in original.groupby(['assay', 'chrom']):
            frame = frame.sort_values(['start', 'end'])
            current = None
            for row in frame.itertuples():
                if current is None or row.start > current['end'] + 1:
                    if current:
                        union.append(current)
                    current = dict(assay=assay, chrom=chrom, start=row.start, end=row.end, genome='mm10')
                else:
                    current['end'] = max(current['end'], row.end)
            if current:
                union.append(current)
        result = pd.DataFrame(union)
        result['peak_id'] = result.assay + ':' + result.chrom + ':' + result.start.astype(str) + '-' + result.end.astype(str)
        result['coordinate_system'] = 'closed1'
        result.to_csv(out / 'mm10_peak_unions.csv', index=False)
    main()
    return locals()

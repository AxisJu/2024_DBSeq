def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import bisect
    import gzip
    import re
    from collections import defaultdict
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def gene_tss(path):
        rows = defaultdict(list)
        op = gzip.open if str(path).endswith('.gz') else open
        with op(path, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                f = line.rstrip().split('\t')
                if len(f) != 9 or f[2] != 'gene':
                    continue
                attrs = dict(re.findall('(\\w+) "([^"]+)"', f[8]))
                tss = int(f[3] if f[6] == '+' else f[4])
                rows[f[0]].append((tss, attrs['gene_id'], attrs.get('gene_name', attrs['gene_id']), f[6]))
        return {chrom: sorted(set(records)) for chrom, records in rows.items()}

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260921-CUT&TAG')
        import pandas as pd
        peaks = inputs['peaks_4'][str(input_path(data, a.peaks))]
        if not peaks.genome.eq('mm10').all() or not peaks.coordinate_system.eq('closed1').all():
            raise ValueError('All peaks must have mm10 closed1 coordinates')
        tss = gene_tss(input_path(data, a.gtf))
        if not tss:
            raise ValueError('No GTF gene features')
        assays = sorted(set(peaks.assay) - {a.assay})
        other = {(assay, chrom): frame.sort_values('start') for (assay, chrom), frame in peaks.groupby(['assay', 'chrom'])}
        rows = []
        for peak in peaks[peaks.assay == a.assay].itertuples():
            records = tss.get(peak.chrom, [])
            if not records:
                raise ValueError('GTF lacks contig: ' + peak.chrom)
            midpoint = (peak.start + peak.end) / 2
            distances = [abs(r[0] - midpoint) for r in records]
            minimum = min(distances)
            nearest = [r for r, d in zip(records, distances) if d == minimum]
            for t, gene_id, gene, strand in nearest:
                row = dict(peak_id=peak.peak_id, chrom=peak.chrom, start=peak.start, end=peak.end, gene_id=gene_id, gene=gene, tss=t, gene_strand=strand, signed_TSS_distance_bp=(midpoint - t) * (1 if strand == '+' else -1))
                for assay in assays:
                    table = other.get((assay, peak.chrom))
                    row['overlap_' + assay] = bool(table is not None and ((table.start <= peak.end) & (table.end >= peak.start)).any())
                rows.append(row)
        annotated = pd.DataFrame(rows)
        if a.trgs:
            trgs = inputs['trgs_5'][str(input_path(data, a.trgs))]
            if trgs.gene.duplicated().any():
                raise ValueError('TRG table must have one row per tested gene')
            annotated = annotated.merge(trgs, on='gene', how='left', validate='many_to_one')
        annotated.to_csv(out / 'peak_TSS_and_interval_overlap.csv', index=False)
        columns = ['overlap_' + assay for assay in assays]
        if columns:
            annotated.groupby(['gene_id', 'gene'])[columns].any().reset_index().to_csv(out / 'gene_interval_overlap.csv', index=False)
    main()
    return locals()

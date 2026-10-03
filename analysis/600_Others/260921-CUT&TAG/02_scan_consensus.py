def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import re
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']
    Fasta = inputs['Fasta_3']
    PATTERNS = {'AP1_TRE': 'TGA[CG]TCA', 'CRE': 'TGACGTCA', 'AP1_variant': 'TGA[CG]TGA', 'GRE_full': 'AGAACA[ATGC]{3}TGTTCT|GGTACA[ATGC]{3}TGTTCT|AGAACA[ATGC]{3}TGTACC|AGAACA[ATGC]{3}TGTCCT', 'GRE_half': 'AGAACA|TGTTCT|GGTACA|TGTACC'}

    def scan(sequence, start1=1):
        reverse = sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
        found = {}
        for name, pattern in PATTERNS.items():
            for strand, seq in [('+', sequence), ('-', reverse)]:
                for match in re.finditer('(?=(' + pattern + '))', seq):
                    a, b = (match.start(1), match.end(1))
                    if strand == '-':
                        a, b = (len(seq) - b, len(seq) - a)
                    key = (name, start1 + a, start1 + b - 1)
                    found.setdefault(key, set()).add(strand)
        return [dict(motif=k[0], start=k[1], end=k[2], strand=''.join(sorted(v))) for k, v in sorted(found.items())]

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_4']
        data, out = roots('600_Others/260921-CUT&TAG')
        import pandas as pd
        peaks = inputs['peaks_5'][str(input_path(data, a.peaks))]
        if not peaks.genome.eq('mm10').all() or not peaks.coordinate_system.eq('closed1').all():
            raise ValueError('Expected mm10 one-based closed peaks')
        if 'assay' in peaks:
            peaks = peaks[peaks.assay == a.assay]
        genome = Fasta(inputs['genome_sequences'])
        hits, rows = ([], [])
        try:
            for peak in peaks.itertuples():
                sequence = genome.fetch(peak.chrom, int(peak.start), int(peak.end))
                found = scan(sequence, int(peak.start))
                row = dict(peak_id=peak.peak_id, chrom=peak.chrom, start=peak.start, end=peak.end, genome='mm10', coordinate_system='closed1')
                for name in PATTERNS:
                    subset = [h for h in found if h['motif'] == name]
                    row['n_' + name] = len(subset)
                for hit in found:
                    hit.update(peak_id=peak.peak_id, chrom=peak.chrom, signed_center_distance_bp=(hit['start'] + hit['end'] - peak.start - peak.end) / 2)
                    hits.append(hit)
                rows.append(row)
        finally:
            genome.close()
        pd.DataFrame(hits).to_csv(out / 'consensus_hits.csv', index=False)
        pd.DataFrame(rows).to_csv(out / 'peak_consensus_counts.csv', index=False)
    main()
    return locals()

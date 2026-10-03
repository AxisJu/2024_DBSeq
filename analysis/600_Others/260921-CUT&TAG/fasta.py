def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    class Fasta:

        def __init__(self, sequences):
            self.sequences = sequences

        def fetch(self, chrom, start, end):
            sequence = self.sequences[chrom]
            length = len(sequence)
            if not 1 <= start <= end <= length:
                raise ValueError('Interval outside FASTA contig')
            return sequence[start - 1:end].upper()

        def close(self):
            return None
    return locals()

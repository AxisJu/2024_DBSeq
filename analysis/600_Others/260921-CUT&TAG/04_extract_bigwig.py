def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260921-CUT&TAG')
        import pandas as pd
        import numpy as np
        tracks = inputs['tracks_4'][str(input_path(data, a.tracks))]
        windows = inputs['windows_5'][str(input_path(data, a.windows))]
        if not tracks.genome.eq('mm10').all() or a.bin_bp < 1:
            raise ValueError('mm10 tracks and a positive bin width are required')
        rows = []
        for track in tracks.itertuples():
            signal = inputs['signal_by_sample'][track.sample_id]
            for window in windows.itertuples():
                if window.chrom not in signal or not 1 <= window.start <= window.end <= len(signal[window.chrom]):
                    raise ValueError('Window is outside the BigWig reference')
                for start0 in range(int(window.start) - 1, int(window.end), a.bin_bp):
                    end0 = min(start0 + a.bin_bp, int(window.end))
                    values = np.asarray(signal[window.chrom][start0:end0], dtype=float)
                    covered = np.isfinite(values)
                    rows.append(dict(sample_id=track.sample_id, window_id=window.window_id, chrom=window.chrom, start=start0 + 1, end=end0, n_covered_bp=int(covered.sum()), mean_signal=float(values[covered].mean()) if covered.any() else np.nan, signal_unit=track.signal_unit, coordinate_system='closed1'))
        pd.DataFrame(rows).to_csv(out / 'bigwig_binned_signal.csv', index=False)
    main()
    return locals()

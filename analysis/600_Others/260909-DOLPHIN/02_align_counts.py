def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import csv
    import subprocess
    from pathlib import Path
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_3']
        data, out = roots('600_Others/260909-DOLPHIN')
        for name in ['03_exon_star', '04_exon_gene_cnt', '05_exon_junct_cnt']:
            (out / name).mkdir(exist_ok=True)
        handle = inputs['handle_4'][str(input_path(data, a.bam_manifest))]
        rows = list(csv.DictReader(handle))
        if len({r['CB'] for r in rows}) != len(rows):
            raise ValueError('Duplicate cell BAM entries')
        junction_manifest = []
        for row in rows:
            cell = row['CB']
            if Path(cell).name != cell or '.aggr' in cell:
                raise ValueError('Expected original cell IDs')
            bam = input_path(data, row['bam'])
            cell_out = out / '03_exon_star' / cell
            cell_out.mkdir(exist_ok=True)
            prefix = str(cell_out / (cell + '.'))
            inputs['data_5']
            aligned = prefix + 'Aligned.sortedByCoord.out.bam'
            for exon, dirname, suffix in [(False, '04_exon_gene_cnt', '.exongene.count.txt'), (True, '05_exon_junct_cnt', '.exon.count.txt')]:
                command = ['featureCounts', '-T', str(a.threads), '-t', 'exon', '-O', '-M']
                if exon:
                    command += ['-f', '-J']
                command += ['-a', str(input_path(data, a.exon_gtf)), '-o', str(out / dirname / (cell + suffix)), aligned]
                inputs['data_6']
            junction_manifest.append({'sample_id': cell, 'sj_path': prefix + 'SJ.out.tab', 'bam': aligned})
        handle = inputs['handle_7'][str(out / 'raw_junction_manifest.csv')]
        writer = csv.DictWriter(handle, fieldnames=['sample_id', 'sj_path', 'bam'])
        writer.writeheader()
        writer.writerows(junction_manifest)
    main()
    return locals()

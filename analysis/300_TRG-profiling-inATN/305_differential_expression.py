def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']

    def main():
        p = parser(__doc__)
        args, extra = inputs['args_extra_3']
        Paths(args)
        import subprocess
        folder = Path(__file__).resolve().parent.parent / '200_TRG-identification'
        mapping = {'edger': '201_edger.R', 'memento': '201_memento.py', 'cossim': '204_cossim.R', 'identify': '205_identify_trg.R'}
        script = folder / mapping[args.method]
        command = [sys.executable if script.suffix == '.py' else args.rscript, str(script), '--data-root', args.data_root, '--output-root', args.output_root, '--scope', 'atn', '--threads', str(args.threads), '--seed', str(args.seed)]
        if script.suffix == '.R':
            command += ['--input', 'data/snRNAseq_mouse/processed/matrix/running_ATN_250704.qs']
        if args.force:
            command.append('--force')
        inputs['data_4']
    main()
    return locals()

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    from pathlib import Path
    import os

    def roots(section):
        data = os.environ.get('DBSEQ_DATA_ROOT')
        output = os.environ.get('DBSEQ_OUTPUT_ROOT')
        if not data or not output:
            raise ValueError('Set DBSEQ_DATA_ROOT and DBSEQ_OUTPUT_ROOT')
        data, output = (Path(data).resolve(), Path(output).resolve())
        public = next((p for p in Path(__file__).resolve().parents if p.name == 'analysis')).parent
        if output == public or public in output.parents:
            raise ValueError('Outputs must be outside the publication code package')
        output = output / section
        output.mkdir(parents=True, exist_ok=True)
        return (data, output)

    def input_path(root, value):
        path = Path(value).expanduser()
        return path.resolve() if path.is_absolute() else (root / path).resolve()
    return locals()

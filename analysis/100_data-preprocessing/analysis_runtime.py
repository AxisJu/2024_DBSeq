def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import os
    import random
    from pathlib import Path
    PACKAGE_ROOT = Path(__file__).resolve().parents[2]

    def parser(description):
        p = argparse.ArgumentParser(description=description)
        return p

    class Paths:

        def __init__(self, args):
            self.data = Path(args.data_root).expanduser().resolve()
            self.out = Path(args.output_root).expanduser().resolve()
            self.external = Path(args.external_root).expanduser().resolve() if args.external_root else self.data / 'external'
            if self.out == PACKAGE_ROOT or PACKAGE_ROOT in self.out.parents:
                raise ValueError('Output root must be outside the publication code package.')
            self.force = args.force
            random.seed(args.seed)
            for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                os.environ[key] = str(args.threads)

        def input(self, relative):
            value = Path(relative)
            if value.is_absolute():
                return value
            if value.parts and value.parts[0] == 'external':
                return self.external.joinpath(*value.parts[1:])
            generated = self.out / value
            return generated if inputs['data_1'] else self.data / value

        def output(self, relative):
            value = (self.out / relative).resolve()
            if value != self.out and self.out not in value.parents:
                raise ValueError('Output paths must remain inside output root.')
            if value == PACKAGE_ROOT or PACKAGE_ROOT in value.parents:
                raise ValueError('Refusing to write inside the code package.')
            value.parent.mkdir(parents=True, exist_ok=True)
            return value
    return locals()

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    from pathlib import Path
    import importlib
    import importlib.util
    import sys

    def load_adaptation(name, namespace='graph_generation'):
        parent = importlib.import_module('DOLPHIN.' + namespace)
        key = 'DOLPHIN.' + namespace + '.' + name
        path = Path(__file__).resolve().parent / 'adaptations' / (name + '.py')
        spec = inputs['spec_1']
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        inputs['data_2']
        setattr(parent, name, module)
        return module

    def patch_original(name, replacements):
        import inspect
        module = importlib.import_module(name)
        source = inputs['source_3']
        for old, new in replacements:
            if old not in source:
                raise RuntimeError('Unsupported upstream version: ' + name)
            source = source.replace(old, new)
        exec(compile(source, module.__file__, 'exec'), module.__dict__)
        return module
    return locals()

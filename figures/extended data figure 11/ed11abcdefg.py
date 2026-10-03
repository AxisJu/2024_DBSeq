def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    import os
    from pathlib import Path

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_1']
        if not a.panels or not set(a.panels) <= set('abcdefg'):
            p.error('--panels must contain only a-g')
        Context = inputs['Context_2']
        context = Context(a)
        if set(a.panels) & set('defg'):
            inputs['data_3']
        if 'a' in a.panels or set(a.panels) & set('bc'):
            panel_a = inputs['panel_a_4']
            panel_b = inputs['panel_b_5']
            panel_c = inputs['panel_c_6']
            for key, callback in [('a', panel_a), ('b', panel_b), ('c', panel_c)]:
                if key in a.panels:
                    callback(context)
        if set(a.panels) & set('defg'):
            panel_d = inputs['panel_d_7']
            panels_ef = inputs['panels_ef_8']
            panel_g = inputs['panel_g_9']
            if 'd' in a.panels:
                panel_d(context)
            if set(a.panels) & set('ef'):
                panels_ef(context)
            if 'g' in a.panels:
                panel_g(context)
    main()
    return locals()

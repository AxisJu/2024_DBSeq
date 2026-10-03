def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']
    load_adaptation = inputs['load_adaptation_3']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_4']
        data, out = roots('600_Others/260909-DOLPHIN')
        import numpy as np
        import random
        import torch
        import pyro
        inputs['data_5'][str('model')]
        training = inputs['training_6'][str('train')]
        random.seed(a.seed)
        np.random.seed(a.seed)
        torch.manual_seed(a.seed)
        pyro.set_rng_seed(a.seed)
        params = dict(gat_channel=[2], nhead=1, gat_dropout=0.1, concat=False, list_gra_enc_hid=[128], gra_p_dropout=0.2, z_dim=30, list_fea_dec_hid=[128], list_adj_dec_hid=[128], lr=0.001, batch=a.batch_size, epochs=a.epochs, kl_beta=0.7, fea_lambda=0.5, adj_lambda=0.5)
        target = out / 'result'
        target.mkdir(exist_ok=True)
        training.run_train('streaming', str(input_path(data, a.features)), str(target), params, torch.device(a.device))
    main()
    return locals()

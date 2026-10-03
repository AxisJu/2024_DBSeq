def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']
    model_arguments = inputs['model_arguments_3']

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_4']
        data, out = roots('600_Others/202607-Squidiff')
        import anndata as ad
        import numpy as np
        import torch
        import random
        path = input_path(data, a.training) if a.training else out / 'training_hvg.h5ad'
        x = inputs['x_5'][str(path)]
        n = x.n_vars
        random.seed(a.seed)
        np.random.seed(a.seed)
        torch.manual_seed(a.seed)
        args = model_arguments(n)
        from Squidiff.train_squidiff import run_training
        from Squidiff.train_util import TrainLoop
        import inspect
        source = inputs['source_6']
        if 'self.resume_checkpoint' not in source or 'os.path.join' not in source:
            raise RuntimeError('This wrapper requires the archived Squidiff checkpoint-directory API')
        model_dir = out / 'models_hvg2000'
        model_dir.mkdir(exist_ok=True)
        args.update(data_path=str(path), batch_size=a.batch_size, lr=0.0001, weight_decay=0.0, ema_rate='0.9999', log_interval=500, save_interval=10000, logger_path=str(model_dir / 'logger'), resume_checkpoint=str(model_dir), use_fp16=False, fp16_scale_growth=0.001, use_ddim=True, schedule_sampler='uniform', microbatch=-1, lr_anneal_steps=a.steps, num_channels=128, dropout=0.0)
        import pandas as pd
        pd.DataFrame({'parameter': list(args), 'value': [str(v) for v in args.values()]}).to_csv(out / 'training_parameters.csv', index=False)
        losses = run_training(args)
        pd.DataFrame({'training_loss': losses}).to_csv(out / 'training_losses.csv', index=False)
    main()
    return locals()

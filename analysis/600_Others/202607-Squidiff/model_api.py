def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import os
    import sys

    def external_source():
        source = os.environ.get('DBSEQ_SQUIDIFF_SOURCE')
        if source:
            from pathlib import Path
            source = Path(source).resolve()
            if not inputs['data_1']:
                raise ValueError('DBSEQ_SQUIDIFF_SOURCE must contain the Squidiff package')

    def model_arguments(n_genes):
        external_source()
        from Squidiff.script_util import model_and_diffusion_defaults
        args = model_and_diffusion_defaults()
        args.update(gene_size=n_genes, output_dim=n_genes, num_layers=3, use_encoder=True, diffusion_steps=1000, noise_schedule='linear', use_drug_structure=False, drug_dimension=1024, comb_num=1)
        return args

    def load_model(checkpoint, n_genes, device):
        import torch
        args = model_arguments(n_genes)
        from Squidiff.script_util import create_model_and_diffusion, args_to_dict, model_and_diffusion_defaults
        model, diffusion = create_model_and_diffusion(**args_to_dict(args, model_and_diffusion_defaults().keys()))
        model.load_state_dict(inputs['data_2'][str(checkpoint)])
        return (model.to(device).eval(), diffusion)

    def encode(model, matrix, device, batch_size=256):
        import torch
        import numpy as np
        from scipy import sparse
        chunks = []
        for start in range(0, matrix.shape[0], batch_size):
            block = matrix[start:start + batch_size]
            block = block.toarray() if sparse.issparse(block) else np.asarray(block)
            with torch.no_grad():
                z = model.encoder(torch.as_tensor(block, dtype=torch.float32, device=device), label=None, drug_dose=None)
            chunks.append(z.detach())
        return torch.cat(chunks)

    def attribution(model, mean_expression, direction, device):
        import torch
        x = torch.tensor(mean_expression, dtype=torch.float32, device=device).reshape(1, -1).requires_grad_(True)
        z = model.encoder(x, label=None, drug_dose=None)
        gradient, = torch.autograd.grad((z * direction).sum(), x)
        return gradient.detach().cpu().numpy()[0]

    def paired_decode(model, diffusion, z, direction, scale, noise):
        import torch

        def decode(code):
            kwargs = dict(z_mod=code.reshape(1, -1).expand(noise.shape[0], -1), x_start=None, drug_dose=None, group=None)
            with torch.no_grad():
                return diffusion.ddim_sample_loop(model, shape=noise.shape, noise=noise.clone(), model_kwargs=kwargs)
        return (decode(z).mean(0).cpu().numpy(), decode(z + scale * direction).mean(0).cpu().numpy())
    return locals()

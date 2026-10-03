def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import argparse
    roots = inputs['roots_1']
    input_path = inputs['input_path_2']
    load_model = inputs['load_model_3']
    encode = inputs['encode_4']
    attribution = inputs['attribution_5']
    paired_decode = inputs['paired_decode_6']

    def sensitivity_scores(attr, effect):
        import numpy as np
        keep = np.isfinite(attr) & np.isfinite(effect) & (abs(attr) > 1e-12) & (abs(effect) > 1e-12)
        x, y = (attr[keep], effect[keep])
        if not len(x):
            return {'reversal_fraction': np.nan, 'weighted_reversal_score': np.nan}
        return {'reversal_fraction': np.mean(x * y < 0), 'weighted_reversal_score': -np.sum(x * y) / np.sum(abs(x * y))}

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        a = inputs['a_7']
        data, out = roots('600_Others/202607-Squidiff')
        import anndata as ad
        import numpy as np
        import pandas as pd
        import torch
        torch.manual_seed(a.seed)
        rng = np.random.default_rng(a.seed)
        features = inputs['features_8'][str(input_path(data, a.features))]
        genes = features.var_names.copy()
        donors = inputs['donors_9'][str(input_path(data, a.donor_map))].set_index('Cell_ID') if a.donor_map else None
        if donors is not None and donors.index.duplicated().any():
            raise ValueError('Duplicate cells in donor map')

        def read(value):
            x = inputs['x_10'][str(input_path(data, value))]
            if not x.var_names.is_unique or not x.obs_names.is_unique or (not set(genes) <= set(x.var_names)):
                raise ValueError('Missing checkpoint genes or duplicated IDs')
            x = x[:, genes].copy()
            donor = donors.reindex(x.obs_names).donor_id if donors is not None else x.obs['donor_id']
            if donor.isna().any() or donor.isin(['unassigned', 'doublet', '']).any():
                raise ValueError('All input cells require a resolved biological donor')
            x.obs['donor_id'] = donor.to_numpy()
            return x
        train = inputs['train_11'][str(a.training)]
        regional = inputs['regional_12'][str(a.regions)]
        eval_data = inputs['eval_data_13'][str(a.evaluation)] if a.evaluation else None
        if eval_data is not None and set(eval_data.obs.donor_id) & set(train.obs.donor_id):
            raise ValueError('Evaluation and training donors overlap')
        model, diffusion = inputs['model_diffusion_14'][str(input_path(data, a.checkpoint))]
        z = encode(model, train.X, a.device)
        masks = {group: np.flatnonzero(train.obs.Group_L2.to_numpy() == group) for group in ['DBS_I', 'Sham_I']}
        if min(map(len, masks.values())) == 0:
            raise ValueError('Training requires DBS_I and Sham_I')
        direction = z[masks['DBS_I']].mean(0) - z[masks['Sham_I']].mean(0)
        if direction.norm() == 0:
            raise ValueError('Zero latent perturbation')
        unit = direction.detach() / direction.norm()
        de = inputs['de_15'][str(input_path(data, a.reference_de))].set_index('X')
        if de.index.duplicated().any():
            raise ValueError('Duplicate genes in differential-expression reference')
        reference = de.reindex(genes)
        significant = (reference.FDR < 0.05).to_numpy()
        effect = reference.logFC.to_numpy()
        groups = {'ATN': train[masks['Sham_I']].copy()}
        for region in sorted(regional.obs.Region.unique()):
            groups[str(region)] = regional[regional.obs.Region == region].copy()
        observations, scores, boot = ([], [], [])
        for region, x in groups.items():
            ids = sorted(x.obs.donor_id.unique())
            donor_means = np.stack([np.asarray(x[x.obs.donor_id == donor].X.mean(axis=0)).ravel() for donor in ids])
            attr = attribution(model, donor_means.mean(axis=0), unit, a.device)
            frame = pd.DataFrame({'region': region, 'gene': genes, 'encoder_gradient': attr, 'reference_log2FC': effect, 'reference_FDR': reference.FDR.to_numpy()})
            if a.decode:
                zr = encode(model, x.X, a.device).mean(0)
                noise = torch.randn(a.draws, len(genes), device=a.device)
                sham, virtual = paired_decode(model, diffusion, zr, direction, a.scale, noise)
                frame['decoded_baseline'] = sham
                frame['decoded_virtual'] = virtual
                frame['decoded_difference'] = virtual - sham
            observations.append(frame)
            scores.append(dict(region=region, n_donors=len(ids), **sensitivity_scores(attr[significant], effect[significant])))
            if len(ids) >= 2:
                for b in range(a.bootstrap):
                    sampled = donor_means[rng.integers(0, len(ids), len(ids))].mean(axis=0)
                    draw = attribution(model, sampled, unit, a.device)
                    boot.append(dict(region=region, replicate=b, **sensitivity_scores(draw[significant], effect[significant])))
        pd.concat(observations).to_csv(out / 'gene_sensitivity_and_predictions.csv', index=False)
        pd.DataFrame(scores).to_csv(out / 'regional_sensitivity_scores.csv', index=False)
        boot = pd.DataFrame(boot)
        boot.to_csv(out / 'donor_bootstrap_scores.csv', index=False)
        if len(boot):
            boot.groupby('region')[['reversal_fraction', 'weighted_reversal_score']].quantile([0.025, 0.975]).to_csv(out / 'donor_bootstrap_intervals.csv')
        if eval_data is not None:
            ez = encode(model, eval_data.X, a.device)
            rows = []
            for donor in sorted(eval_data.obs.donor_id.unique()):
                sel = (eval_data.obs.donor_id == donor).to_numpy()
                row = dict(donor_id=donor, n_cells=int(sel.sum()), evaluation='independent_donor_latent_projection')
                row['projection'] = (ez[sel] * unit).sum(1).mean().item()
                rows.append(row)
            pd.DataFrame(rows).to_csv(out / 'independent_donor_projections.csv', index=False)
        pd.DataFrame([dict(scale=a.scale, seed=a.seed, n_training_cells=train.n_obs, evaluation='calibration_and_sensitivity' if eval_data is None else 'independent_donor_projection', score_unit='encoder_gradient', perturbation_unit='model_expression_difference', legacy_scale_selection='same_cohort_calibration')]).to_csv(out / 'analysis_parameters.csv', index=False)
    main()
    return locals()

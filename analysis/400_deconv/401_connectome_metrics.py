def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    Paths = inputs['Paths_1']
    parser = inputs['parser_2']
    import numpy as np

    def shortest_paths(lengths):
        from scipy.sparse.csgraph import shortest_path
        return shortest_path(lengths, directed=True, unweighted=False)

    def diffusion_efficiency(weights):
        from scipy.sparse.csgraph import shortest_path
        weights = np.asarray(weights, dtype=float)
        size = len(weights)
        totals = weights.sum(axis=1, keepdims=True)
        transition = np.divide(weights, totals, out=np.zeros_like(weights), where=totals > 0)
        output = np.zeros_like(weights)
        for target in range(size):
            adjacency = weights > 0
            adjacency[target, :] = False
            reach = np.isfinite(shortest_path(adjacency.astype(float), directed=True, unweighted=True))
            can_reach = reach[:, target]
            bad = ~can_reach
            fails = reach[:, bad].any(axis=1) if bad.any() else np.zeros(size, dtype=bool)
            valid = np.flatnonzero(can_reach & ~fails & (np.arange(size) != target))
            if valid.size:
                q = transition[np.ix_(valid, valid)]
                hitting = np.linalg.solve(np.eye(len(valid)) - q, np.ones(len(valid)))
                output[valid, target] = 1 / hitting
        return output

    def to_length(W, eps=1e-12):
        W = np.array(W, dtype=float)
        L = np.where(W > eps, 1.0 / np.maximum(W, eps), np.inf)
        np.fill_diagonal(L, 0.0)
        return L

    def shortest_path_efficiency(W):
        L = to_length(W)
        D = shortest_paths(L)
        SPE = np.zeros_like(D, dtype=float)
        mask = np.isfinite(D) & (D > 0)
        SPE[mask] = 1.0 / D[mask]
        np.fill_diagonal(SPE, 0.0)
        return SPE

    def navigation_efficiency_from_ED(W, ED, max_hops=None):
        W = np.asarray(W, dtype=float)
        L = to_length(W)
        ED = np.asarray(ED, dtype=float)
        N = W.shape[0]
        if max_hops is None:
            max_hops = 2 * N
        NE = np.zeros((N, N), dtype=float)
        neighbors = [np.where(W[u, :] > 0)[0] for u in range(N)]
        for i in range(N):
            for j in range(N):
                if i == j:
                    continue
                u = i
                visited = {u}
                cost = 0.0
                success = False
                for _ in range(max_hops):
                    nbrs = neighbors[u]
                    if nbrs.size == 0:
                        break
                    v = nbrs[np.argmin(ED[nbrs, j])]
                    if v in visited or not np.isfinite(L[u, v]):
                        break
                    cost += L[u, v]
                    u = v
                    if u == j:
                        success = True
                        break
                    visited.add(u)
                if success and cost > 0:
                    NE[i, j] = 1.0 / cost
        return NE

    def search_information(W, tol=1e-09):
        W = np.array(W, float)
        N = W.shape[0]
        L = to_length(W)
        D = shortest_paths(L)
        rowSums = W.sum(axis=1, keepdims=True)
        T = np.divide(W, rowSums, out=np.zeros_like(W), where=rowSums > 0)
        SI = np.full((N, N), np.inf, float)
        for i in range(N):
            for j in range(N):
                if i == j or not np.isfinite(D[i, j]):
                    continue
                u = i
                p = 1.0
                ok = True
                for _ in range(N):
                    nbrs = np.where(W[u, :] > 0)[0]
                    eligible = [k for k in nbrs if np.isfinite(L[u, k]) and np.isfinite(D[k, j]) and (abs(L[u, k] + D[k, j] - D[u, j]) < tol)]
                    if len(eligible) == 0:
                        ok = False
                        break
                    k = max(eligible, key=lambda v: T[u, v])
                    p *= max(T[u, k], 1e-300)
                    u = k
                    if u == j:
                        break
                if ok and u == j:
                    SI[i, j] = -np.log2(p)
        return SI

    def communicability(W, normalize=True):
        from scipy.linalg import expm
        W = np.array(W, float)
        if normalize:
            s_out = W.sum(axis=1)
            s_in = W.sum(axis=0)
            s_out_safe = np.where(s_out > 0, s_out, 1.0)
            s_in_safe = np.where(s_in > 0, s_in, 1.0)
            S = np.sqrt(np.outer(s_out_safe, s_in_safe))
            Wp = np.where(W > 0, W / S, 0.0)
        else:
            Wp = W
        CMY = expm(Wp)
        np.fill_diagonal(CMY, 0.0)
        return CMY

    def main():
        p = parser(__doc__)
        args = inputs['args_3']
        paths = Paths(args)
        import pandas as pd
        weights = inputs['weights_4'][str(paths.input(args.weights))]
        if not weights.index.is_unique or not weights.columns.is_unique or set(weights.index) != set(weights.columns):
            raise ValueError('Weights must be square with unique matching row and column labels')
        weights = weights.loc[weights.index, weights.index]
        aligned = []
        for filename in (args.pvalues, args.distances):
            table = inputs['table_5'][str(paths.input(filename))].reindex(index=weights.index, columns=weights.columns)
            if table.isna().any().any():
                raise ValueError('Missing values or unmatched region labels')
            aligned.append(table.to_numpy(float))
        pvalues, distances = aligned
        raw = weights.to_numpy(float)
        if not np.isfinite(raw).all() or (raw < 0).any() or (distances < 0).any():
            raise ValueError('Connectivity and distances must be finite and nonnegative')
        raw = np.where(pvalues <= args.p_threshold, raw, 0)
        transformed = np.zeros_like(raw)
        nonzero = raw > 0
        if nonzero.any():
            logged = np.log10(raw[nonzero])
            transformed[nonzero] = logged + abs(logged.min()) + 1
        outputs = {'SPE': shortest_path_efficiency(transformed), 'NE': navigation_efficiency_from_ED(transformed, distances), 'DE': diffusion_efficiency(transformed), 'SI': search_information(transformed), 'CMY': communicability(transformed)}
        for name, values in outputs.items():
            pd.DataFrame(values, index=weights.index, columns=weights.columns).to_csv(paths.output(f'data/derivatives/dbseq/allen_mousebrainconnectome_{name}_v2.csv'))
    main()
    return locals()

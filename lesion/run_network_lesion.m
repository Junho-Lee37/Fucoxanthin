function run_network_lesion(cfg)
%RUN_NETWORK_LESION

S = load(cfg.fitMatFile);

need = {'alpha','delta','gamma','W','s','K','isTF','hill_n'};
for k = 1:numel(need)
    if ~isfield(S, need{k})
        error('fit_result.mat missing field: %s', need{k});
    end
end

alpha  = S.alpha(:);
delta  = S.delta(:);
gamma  = S.gamma(:);
W0     = S.W;
s      = S.s(:);
K      = S.K(:);
isTF   = logical(S.isTF(:));
hill_n = S.hill_n;

G = numel(alpha);

genes = read_gene_list(cfg.genesFile, G);

if isempty(cfg.lesionTargets)
    lesionTargets = find(isTF);
else
    lesionTargets = cfg.lesionTargets(:);
end

if cfg.useNonTFOnly
    idxOut = find(~isTF);
else
    idxOut = (1:G).';
end

theta0 = pack_theta(alpha, delta, gamma, W0, s, K, isTF, hill_n);

%% -------- baseline (ONCE) --------
if cfg.doBaseline
    [t_none0, x_none0] = simulate_condition(theta0, 'None', cfg.tspan);
    [t_fx0,   x_fx0]   = simulate_condition(theta0, 'Fx',   cfg.tspan);

    base = quantify_effect(t_none0, x_none0, t_fx0, x_fx0, ...
                           cfg.aucWindow, idxOut, cfg.nGrid);
    base.genes  = genes;
    base.idxOut = idxOut;

    % store trajectories for multi-panel plots
    base.t_none = t_none0; base.x_none = x_none0;
    base.t_fx   = t_fx0;   base.x_fx   = x_fx0;

    base.gamma = gamma;

    save(fullfile(cfg.outDir,'baseline.mat'), '-struct', 'base');
end

%% -------- lesion loop --------
results = struct();
results.genes         = genes;
results.isTF          = isTF;
results.idxOut        = idxOut;
results.lesionTargets = lesionTargets;

for m = 1:numel(lesionTargets)
    j = lesionTargets(m);

    W = W0;
    W(:, j) = 0;   % outgoing lesion

    theta = pack_theta(alpha, delta, gamma, W, s, K, isTF, hill_n);

    [t_none, x_none] = simulate_condition(theta, 'None', cfg.tspan);
    [t_fx,   x_fx]   = simulate_condition(theta, 'Fx',   cfg.tspan);

    q = quantify_effect(t_none, x_none, t_fx, x_fx, ...
                        cfg.aucWindow, idxOut, cfg.nGrid);

    results.perLesion(m).lesionedIndex = j;
    results.perLesion(m).lesionedGene  = char(genes(j));
    results.perLesion(m).deltaAUC      = q.deltaAUC;
    results.perLesion(m).AUC_none      = q.AUC_none;
    results.perLesion(m).AUC_fx        = q.AUC_fx;

    % store trajectories
    results.perLesion(m).t_none = t_none;
    results.perLesion(m).x_none = x_none;
    results.perLesion(m).t_fx   = t_fx;
    results.perLesion(m).x_fx   = x_fx;

    fprintf('[%d/%d] lesion %s (index %d) done.\n', ...
            m, numel(lesionTargets), char(genes(j)), j);
end

save(fullfile(cfg.outDir,'lesion_results.mat'), 'results');
end


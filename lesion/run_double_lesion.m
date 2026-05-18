function results2 = run_double_lesion(cfg)
% Double TF lesion: compute interaction I_{i,jk}

Sfit = load(cfg.fitMatFile);

alpha  = Sfit.alpha(:);
delta  = Sfit.delta(:);
gamma  = Sfit.gamma(:);
W0     = Sfit.W;
s      = Sfit.s(:);
K      = Sfit.K(:);
isTF   = logical(Sfit.isTF(:));
hill_n = Sfit.hill_n;

G = numel(alpha);
genes = read_gene_list(cfg.genesFile, G);

% outputs (nonTF only, or all)
if cfg.useNonTFOnly
    idxOut = find(~isTF);
else
    idxOut = (1:G).';
end

theta0 = pack_theta(alpha, delta, gamma, W0, s, K, isTF, hill_n);

% --- baseline ΔAUC ---
[t_none0, x_none0] = simulate_condition(theta0,'None',cfg.tspan);
[t_fx0,   x_fx0]   = simulate_condition(theta0,'Fx',  cfg.tspan);
q0 = quantify_effect(t_none0,x_none0,t_fx0,x_fx0,cfg.aucWindow,idxOut,cfg.nGrid);
dAUC_base = q0.deltaAUC(:);     % (nOut x 1)

% --- single lesion ΔΔAUC (reuse if already computed and saved) ---
% If you already have lesion_results.mat, load it and build dAUC_single(:,m)
L = load(fullfile(cfg.outDir,'lesion_results.mat'));
R = L.results;
lesions = R.lesionTargets(:);
nLes = numel(lesions);
nOut = numel(idxOut);

dAUC_single = zeros(nOut,nLes);
for m=1:nLes
    dAUC_single(:,m) = R.perLesion(m).deltaAUC(:) - dAUC_base;
end

% --- double lesion loop ---
pairs = nchoosek(1:nLes,2);
nPairs = size(pairs,1);

results2.genes = genes;
results2.idxOut = idxOut;
results2.lesionTargets = lesions;
results2.pairs = pairs;

results2.I_gene_pair = zeros(nOut, nPairs);   % I_{i,jk}
results2.I_mean = zeros(nPairs,1);            % mean_i I_{i,jk}
results2.I_meanAbs = zeros(nPairs,1);         % mean_i |I_{i,jk}|

for p=1:nPairs
    a = pairs(p,1);
    b = pairs(p,2);
    j = lesions(a);
    k = lesions(b);

    W = W0;
    W(:, [j k]) = 0;

    theta = pack_theta(alpha, delta, gamma, W, s, K, isTF, hill_n);

    [t_none, x_none] = simulate_condition(theta,'None',cfg.tspan);
    [t_fx,   x_fx]   = simulate_condition(theta,'Fx',  cfg.tspan);

    q = quantify_effect(t_none,x_none,t_fx,x_fx,cfg.aucWindow,idxOut,cfg.nGrid);
    dAUC_double = q.deltaAUC(:);  % (nOut x 1)

    dd_double = dAUC_double - dAUC_base;
    dd_add    = dAUC_single(:,a) + dAUC_single(:,b);

    I = dd_double - dd_add;

    results2.I_gene_pair(:,p) = I;
    results2.I_mean(p)    = mean(I);
    results2.I_meanAbs(p) = mean(abs(I));

    fprintf('[%d/%d] double lesion (%s,%s) done.\n', p,nPairs,genes{j},genes{k});
end

save(fullfile(cfg.outDir,'double_lesion_results.mat'), 'results2');
end

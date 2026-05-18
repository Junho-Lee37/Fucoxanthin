%% ===================== Figure 7 (3x3, FINAL): TF threshold-gated dynamics (gene-order panels) =====================
% Requires (same folder as this script):
%   - fit_result.mat
%   - genes_used.txt
%
% Output:
%   - Fig7_outputs/Fig7_3x3_geneOrder_final.png

clear; clc; close all;

%% ---------- Files ----------
thisDir   = fileparts(mfilename('fullpath'));
fit_file  = fullfile(thisDir, 'fit_result.mat');
gene_file = fullfile(thisDir, 'genes_used.txt');

assert(exist(fit_file,'file')==2,  'Cannot find fit_result.mat in: %s', thisDir);
assert(exist(gene_file,'file')==2, 'Cannot find genes_used.txt in: %s', thisDir);

%% ---------- Style (requested) ----------
FS_LABEL = 25;   % axis labels
FS_TICK  = 20;   % tick numbers
LW_LINE  = 3;    % curves/lines
LW_AX    = 1.5;  % axes box
FS_SGT   = 22;   % sgtitle font

XTICKS = [0 40 80 120 160];

% Colors
c_state = [0 0 0];          % TF state curve
c_thr   = [0.35 0.35 0.35]; % threshold line
c_onset = [0.20 0.20 0.20]; % fuco onset line

%% ---------- Time / ODE ----------
t1 = 0:1:80;
t2 = 80:1:160;
ode_opts   = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',2);
ode_x0_cap = 1e3;

%% ---------- Load fitted result ----------
S = load(fit_file);

W      = S.W;
alpha  = S.alpha(:);
delta  = S.delta(:);
gamma  = S.gamma(:);
s      = S.s(:);
K      = S.K(:);
isTF   = logical(S.isTF(:));
tf_idx = S.tf_idx(:);              % gene indices of TFs (order matters!)
hill_n = double(S.hill_n);

G   = numel(alpha);
nTF = numel(tf_idx);

% Gene names in model order
genes = readlines(gene_file);
genes = strtrim(genes);
genes(genes=="") = [];
assert(numel(genes)==G, 'genes_used.txt count (%d) must match model size G (%d).', numel(genes), G);

% Map gene index -> TF threshold index
tf_map = zeros(G,1);
tf_map(tf_idx) = 1:nTF;

%% ---------- Select 7 TFs in GENE ORDER ----------
nShow = min(7, nTF);
tf_show_idx = tf_idx(1:nShow);   % <- gene order 그대로

%% ---------- Initial condition ----------
x0 = min(max(s, 1e-3), ode_x0_cap);

%% ---------- RHS (matches fitting code logic) ----------
rhs = @(t,x,Fxflag) rhs_fitstyle(t,x,Fxflag,G,W,alpha,delta,gamma,isTF,tf_map,K,hill_n);

%% ---------- Piecewise simulation (0–80 OFF, 80–160 ON) ----------
[~, X1] = ode15s(@(t,x) rhs(t,x,0), t1, x0, ode_opts);
x80 = X1(end,:).';
[~, X2] = ode15s(@(t,x) rhs(t,x,1), t2, x80, ode_opts);

t_all = [t1(:); t2(2:end).'];
X_all = [X1;     X2(2:end,:)];

%% ---------- Metrics for summary panels (in the SAME gene order) ----------
out_show   = zeros(nShow,1);  % out-strength Σ_i |W_ij|
prox_after = zeros(nShow,1);  % min_{t>=80} |x/K - 1|
mask_after = (t_all >= 80);

for m = 1:nShow
    j = tf_show_idx(m);
    kpos = tf_map(j);
    out_show(m) = sum(abs(W(:,j)));
    ratio = X_all(:,j) ./ max(K(kpos), 1e-12);
    prox_after(m) = min(abs(ratio(mask_after) - 1));
end

%% ===================== Plot (3x3) =====================
fig = figure('Color','w','Position',[50 50 1550 1500]); % balanced aspect
tl  = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

% ---- Panels 1–7: TF trajectories x_j(t) with K_j + onset ----
for m = 1:nShow
    j = tf_show_idx(m);
    kpos = tf_map(j);

    ax = nexttile(m); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    plot(ax, t_all, X_all(:,j), 'LineWidth', LW_LINE);
    yline(ax, K(kpos), '--', 'LineWidth', LW_LINE, 'Color', c_thr);
    xline(ax, 80,      '-',  'LineWidth', LW_LINE, 'Color', c_onset);

    ylabel(ax, sprintf('%s\nx_j(t)', genes(j)), ...
        'Interpreter','none', 'FontSize', FS_LABEL, 'FontWeight','bold');

    if m >= 7
        xlabel(ax, 'Time (min)', 'FontSize', FS_LABEL, 'FontWeight','bold');
    end

    xlim(ax, [0 160]);
    xticks(ax, XTICKS);
    set(ax, 'FontSize', FS_TICK, 'LineWidth', LW_AX);
end

% ---- Panel 8: out-strength bar (GENE ORDER, aligned with panels 1–7) ----
ax8 = nexttile(8); hold(ax8,'on'); box(ax8,'on'); grid(ax8,'on');
bar(ax8, 1:nShow, out_show, 'LineWidth', 1.0);
xticks(ax8, 1:nShow);
xticklabels(ax8, genes(tf_show_idx));
xtickangle(ax8, 45);
axis([0.5 7.5 0 0.3]);

ylabel(ax8, {'TF out-strength','\Sigma_i |W_{ij}|'}, ...
    'FontSize', FS_LABEL, 'FontWeight','bold');
xlabel(ax8, 'TF', 'FontSize', FS_LABEL, 'FontWeight','bold');

set(ax8, 'FontSize', FS_TICK, 'LineWidth', LW_AX);

% ---- Panel 9: threshold proximity scatter (points labeled; same order) ----
ax9 = nexttile(9); hold(ax9,'on'); box(ax9,'on'); grid(ax9,'on');
scatter(ax9, out_show, prox_after, 150, 'filled');

xlabel(ax9, 'TF out-strength  \Sigma_i |W_{ij}|', ...
    'FontSize', FS_LABEL, 'FontWeight','bold');
ylabel(ax9, {'Threshold proximity', 'min_{t\in[80,160]} |x_j/K_j - 1|'}, ...
    'FontSize', FS_LABEL, 'FontWeight','bold');

set(ax9, 'FontSize', FS_TICK, 'LineWidth', LW_AX);

for m = 1:nShow
    text(ax9, out_show(m), prox_after(m), "  " + genes(tf_show_idx(m)), ...
        'Interpreter','none', 'FontSize', 16, 'FontWeight','bold');
end

% ---- Global title ----
sgtitle(tl, 'Figure 4 | Threshold-gated TF dynamics under fucoxanthin', ...
    'FontSize', FS_SGT, 'FontWeight','bold');

%% ---------- FORCE FONT OVERRIDE (prevents tiledlayout/export resets) ----------
axs = findall(fig, 'Type', 'axes');
for k = 1:numel(axs)
    ax = axs(k);
    ax.FontSize  = FS_TICK;  % tick numbers
    ax.LineWidth = LW_AX;    % box

    if ~isempty(ax.XLabel)
        ax.XLabel.FontSize = FS_LABEL;
        ax.XLabel.FontWeight = 'bold';
    end
    if ~isempty(ax.YLabel)
        ax.YLabel.FontSize = FS_LABEL;
        ax.YLabel.FontWeight = 'bold';
    end

    % Ensure TF panels have the requested xticks (summary panels excluded automatically)
    % (Only apply if axis range matches TF panels roughly)
    xl = ax.XLim;
    if xl(1) <= 0 && xl(2) >= 160 && numel(ax.XTick) > 0
        ax.XTick = XTICKS;
    end
end

%% ---------- Save ----------
outDir = fullfile(thisDir, 'Fig7_outputs');
if ~exist(outDir,'dir'), mkdir(outDir); end

outPng = fullfile(outDir, 'Fig7_3x3_geneOrder_final.png');
exportgraphics(fig, outPng, 'Resolution', 300);
fprintf('✅ Saved: %s\n', outPng);

%% ===================== Local function: RHS (matches fitting code) =====================
function dx = rhs_fitstyle(~, x, Fxflag, G, W, alpha, delta, gamma, isTF, tf_map, K, hill_n)
    dx = alpha - delta.*x;
    if Fxflag ~= 0
        dx = dx + gamma;
    end

    for i = 1:G
        acc = 0;
        for j = 1:G
            if j == i, continue; end
            Wij = W(i,j);
            if Wij == 0, continue; end

            xj = x(j);
            if isTF(j)
                kidx = tf_map(j);
                Kj   = K(kidx);
                xj_pos = max(xj, 0);
                xpow = xj_pos^hill_n;
                den  = (Kj^hill_n) + xpow + 1e-12;
                h    = xpow / den;           % in [0,1]
                acc  = acc + Wij * h;
            else
                acc  = acc + Wij * xj;       % linear for non-TF sources
            end
        end
        dx(i) = dx(i) + acc;
    end
end

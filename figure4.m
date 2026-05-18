%% ===================== Fig.3 (ALL genes, 0–80): RHS matches fitting code =====================
% Matches your fitting code exactly:
%  - RHS: dx = alpha - delta.*x + (Fxflag)*gamma + W*(TF Hill / nonTF linear)
%  - s is NOT in RHS; s is applied ONLY to observed data when plotting (s ∘ data)
%  - Add t=0 data point as s_i using SAME markers (no special marker)
%  - x0 used for plotting: x0 = clamp(s) (consistent with your "x0_s" usage)
%  - Fuco OFF vs ON: Fxflag = 0/1 over 0–80
%  - Plot ALL genes from genes_used.txt (order consistent with fit_result.mat)
%
% Required in current folder:
%   fit_result.mat, genes_used.txt, data2.xlsx

clear; clc; close all;

%% ---- Files ----
fit_file   = 'fit_result.mat';
gene_file  = 'genes_used.txt';
xls_file   = 'data2.xlsx';

%% ---- Time / solver options (match your fitting script) ----
tsamp     = [10 20 40 80];
tgrid     = 0:1:80;
ode_opts  = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',2);

% For data points (include t=0 as data point)
data_t = [0 tsamp];

%% ---- Layout / output ----
nRow = 5; nCol = 4;       % 20 genes per page
perPage = nRow*nCol;

outDir = 'Fig3_all_genes_pages';
if ~exist(outDir,'dir'), mkdir(outDir); end
savePNG = true; pngDPI = 200;

%% ---- x0 cap (from your script) ----
ode_x0_cap = 1e3;

%% ===================== 1) Load fitted parameters =====================
S = load(fit_file);

W      = S.W;
alpha  = S.alpha(:);
delta  = S.delta(:);
gamma  = S.gamma(:);
s      = S.s(:);
K      = S.K(:);                 % length nTF
isTF   = logical(S.isTF(:));     % length G
tf_idx = S.tf_idx(:);            % TF indices (length nTF)
hill_n = double(S.hill_n);

G = numel(alpha);

% Build tf_map: gene index j -> K index (0 if nonTF), EXACTLY like your fitting code
tf_map = zeros(G,1);
tf_map(tf_idx) = 1:numel(tf_idx);
nTF = numel(tf_idx);

%% ===================== 2) Load gene list (model order) =====================
geneNames = readlines(gene_file);
geneNames = strtrim(geneNames);
geneNames(geneNames=="") = [];

assert(numel(geneNames)==G, 'genes_used.txt count must match fit_result.mat gene count.');

%% ===================== 3) Load data2.xlsx and reconstruct Fuco absolute =====================
none_cols = {'IL1b_10m/None.fc','IL1b_20m/None.fc','IL1b_40m/None.fc','IL1b_80m/None.fc'};
fx_cols   = {'IL1b_Fx_10m/IL1b_10m.fc','IL1b_Fx_20m/IL1b_20m.fc','IL1b_Fx_40m/IL1b_40m.fc','IL1b_Fx_80m/IL1b_80m.fc'};

T = readtable(xls_file,'VariableNamingRule','preserve');
vars = T.Properties.VariableNames;

% detect gene column (same logic as yours)
gene_col = '';
cands = {'Gene_Symbol','Gene','Symbol','GeneID','Transcript_ID','MGI'};
for i=1:numel(cands)
    if ismember(cands{i}, vars), gene_col = cands{i}; break; end
end
if isempty(gene_col), gene_col = vars{1}; end

genes_all = string(T.(gene_col));
genes_all = strtrim(genes_all);

assert(all(ismember(none_cols, vars)), 'None columns missing in data2.xlsx');
assert(all(ismember(fx_cols,   vars)), 'Fx ratio columns missing in data2.xlsx');

none_expr_all = table2array(T(:, none_cols));   % [G_all x 4]
fx_ratio_all  = table2array(T(:, fx_cols));     % [G_all x 4]

% Negative handling: -x -> 1/|x|
none_expr_all(none_expr_all<0) = 1./abs(none_expr_all(none_expr_all<0));
fx_ratio_all( fx_ratio_all <0) = 1./abs(fx_ratio_all( fx_ratio_all <0));

fx_expr_all = none_expr_all .* fx_ratio_all;    % Fuco absolute

% Match excel rows to model gene order
[tfGene, loc] = ismember(geneNames, genes_all);

Y_il1b = nan(G,4);
Y_fuc  = nan(G,4);
present = find(tfGene);
rows    = loc(present);

Y_il1b(present,:) = none_expr_all(rows,:);
Y_fuc(present,:)  = fx_expr_all(rows,:);

fprintf('[Fig3] Data matched: present %d / %d genes (missing -> NaN)\n', numel(present), G);

%% ===================== 4) Apply s to observed data (s ∘ data), robustly =====================
% This is EXACTLY what your residual does: data_scaled = data .* s_row
Y_il1b_s = bsxfun(@times, Y_il1b, s);   % [G x 4]
Y_fuc_s  = bsxfun(@times, Y_fuc,  s);   % [G x 4]

%% ===================== 5) ODE simulation (RHS consistent with fitting code) =====================
% Use x0 = clamp(s) (matches your plotting choice "x0_s")
x0 = min(max(s(:), 1e-3), ode_x0_cap);

rhs = @(t,x,Fxflag) rhs_fitstyle(t,x, Fxflag, G, W, alpha, delta, gamma, isTF, tf_map, K, hill_n);

% Solve on the same grid you used (0:1:80) so interpolation matches
[~, XN] = ode15s(@(t,x) rhs(t,x,0), tgrid, x0, ode_opts);  % Fuco OFF
[~, XF] = ode15s(@(t,x) rhs(t,x,1), tgrid, x0, ode_opts);  % Fuco ON

% Sample model at tsamp (same as your code; use pchip)
YN_samp = interp1(tgrid, XN, tsamp, 'pchip');   % [4 x G]
YF_samp = interp1(tgrid, XF, tsamp, 'pchip');   % [4 x G]

%% ===================== 6) Plot ALL genes across pages =====================
nPages = ceil(G / perPage);

for pg = 1:nPages
    figure('Color','w','Position',[60 60 1400 900]);

    startIdx = (pg-1)*perPage + 1;
    endIdx   = min(pg*perPage, G);

    for idx = startIdx:endIdx
        sp = idx - startIdx + 1;
        subplot(nRow, nCol, sp); hold on;

        % --- Model curves (same units as s ∘ data) ---
        plot(tgrid, XN(:,idx), 'LineWidth', 3);   % F=0
        plot(tgrid, XF(:,idx), 'LineWidth', 3);   % F=1

        % --- Data points: [t=0,10,20,40,80] with SAME markers ---
        if all(isfinite(Y_il1b_s(idx,:)))
            y_il1b = [s(idx), Y_il1b_s(idx,:)];     % t=0 is s_i
            plot(data_t, y_il1b, 'o', 'MarkerSize', 10, 'LineWidth', 1.5);
        else
            text(0.02, 0.90, 'IL1b data missing', 'Units','normalized', 'FontSize',10, 'FontAngle','italic');
        end

        if all(isfinite(Y_fuc_s(idx,:)))
            y_fuc = [s(idx), Y_fuc_s(idx,:)];
            plot(data_t, y_fuc,  's', 'MarkerSize', 10, 'LineWidth', 1.5);
        else
            text(0.02, 0.80, 'Fuco data missing', 'Units','normalized', 'FontSize',10, 'FontAngle','italic');
        end

        title(geneNames(idx), 'Interpreter','none', 'FontWeight','normal');
        xlim([0 80]); box on;

        if sp > (perPage - nCol), xlabel('Time (min)'); end
        if mod(sp-1, nCol) == 0, ylabel('Expression (a.u.)'); end

        if sp == 1
            legend({'Model F=0 (IL-1\beta)','Model F=1 (IL-1\beta+Fuco)', ...
                    'Data IL-1\beta (s·data, t0=s)','Data IL-1\beta+Fuco (s·data, t0=s)'}, ...
                   'Location','best');
        end
    end

    sgtitle(sprintf('Figure 3 (all genes, 0–80 min) — page %d/%d', pg, nPages), 'FontWeight','normal');
    % ===== Global font scaling for this figure =====
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    set(findall(gcf,'-property','FontName'),'FontName','Arial')

    if savePNG
        outPng = fullfile(outDir, sprintf('Fig3_page_%02d.png', pg));
        exportgraphics(gcf, outPng, 'Resolution', pngDPI);
    end
end

fprintf('[Fig3] Saved %d page(s) to: %s\n', nPages, outDir);

%% ===================== Local RHS (identical logic to local_rhs_all2all_TF) =====================
function dx = rhs_fitstyle(~, x, Fxflag, G, W, alpha, delta, gamma, isTF, tf_map, K, hill_n)
    % Base
    dx = alpha - delta.*x;
    if Fxflag ~= 0
        dx = dx + gamma;
    end

    % Coupling: source j -> target i (W(i,j))
    for i = 1:G
        acc = 0;
        for j = 1:G
            if j==i, continue; end
            Wij = W(i,j);
            if Wij == 0, continue; end

            xj = x(j);
            if isTF(j)
                kidx = tf_map(j);        % 1..nTF
                Kj   = K(kidx);
                xj_pos = max(xj, 0);
                xpow = xj_pos^hill_n;
                den  = (Kj^hill_n) + xpow + 1e-12;
                h    = xpow / den;       % in [0,1]
                acc  = acc + Wij * h;
            else
                acc  = acc + Wij * xj;   % linear for non-TF sources
            end
        end
        dx(i) = dx(i) + acc;
    end
end


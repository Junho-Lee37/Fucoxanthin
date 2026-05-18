%% =========================================================
% Fig.5 (Downstream non-TF genes): FIXED axes positions (Fig1==Fig2 spacing)
%  - Fig1: 3×3 non-TF #1–9
%  - Fig2: SAME 3×3 geometry
%      tiles 1–4: non-TF #10–13
%      tile 5: ΔAUC bar
%      tile 6: divergence heatmap (80–160) + manual colorbar
%      tiles 7–9: blank
%
% Key: axes are placed by explicit Position -> identical gaps across figures.
%% =========================================================
clear; clc; close all;

%% ---- Global font defaults (Times New Roman) ----
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultLegendFontName','Times New Roman');
set(groot,'defaultColorbarFontName','Times New Roman');

%% ---- Files (assume same folder) ----
fit_file   = 'fit_result.mat';
genes_file = 'genes_used.txt';
assert(exist(fit_file,'file')==2, 'Cannot find %s', fit_file);
assert(exist(genes_file,'file')==2, 'Cannot find %s', genes_file);

S = load(fit_file);
genes = string(readlines(genes_file)); genes = strtrim(genes); genes = genes(genes~="");

W      = S.W;
alpha  = S.alpha(:);
delta  = S.delta(:);
gamma  = S.gamma(:);
s      = S.s(:);
isTF   = logical(S.isTF(:));
K      = S.K(:);
hill_n = S.hill_n;

G = size(W,1);
if numel(genes) ~= G
    warning('genes_used length (%d) != size(W,1) (%d). Adjusting.', numel(genes), G);
    if numel(genes) > G, genes = genes(1:G);
    else, genes(end+1:G) = "GENE_" + string(numel(genes)+1:G);
    end
end

nonTF_idx   = find(~isTF);
genes_nonTF = genes(nonTF_idx);
N = numel(nonTF_idx);
assert(N==13, 'Expected 13 non-TF genes, got %d.', N);

%% ---- Time grids & ODE options ----
tgrid = (0:1:160).';
t1    = (0:1:80).';
t2    = (80:1:160).';
ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',2);
xt_all = [0 40 80 120 160];

%% ---- Initial condition: x0 = s ----
x0 = max(s, 1e-3);

%% ---- RHS (TF=Hill, nonTF=linear) ----
rhs = @(t,x,Fxflag) local_rhs(x, Fxflag, alpha, delta, gamma, W, isTF, K, hill_n);

%% ---- Simulations ----
% SWITCH (0–80 noFx, 80–160 Fx)
[~, X_0to80] = ode15s(@(t,x) rhs(t,x,0), t1, x0, ode_opts);
x80 = X_0to80(end,:).';
[~, X_80to160] = ode15s(@(t,x) rhs(t,x,1), t2, x80, ode_opts);
X_switch = [X_0to80; X_80to160(2:end,:)];  % 161×G

% Counterfactual (noFx, 0–160)
[~, X_noFx] = ode15s(@(t,x) rhs(t,x,0), tgrid, x0, ode_opts); % 161×G

% ΔAUC (80–160) for non-TF
i80 = find(tgrid==80,1);
dAUC = zeros(N,1);
for k=1:N
    g = nonTF_idx(k);
    dAUC(k) = trapz(tgrid(i80:end), X_switch(i80:end,g) - X_noFx(i80:end,g));
end

% Heatmap data (80–160)
t_post = tgrid(i80:end);
D = X_switch(i80:end, nonTF_idx) - X_noFx(i80:end, nonTF_idx); % 81×13

%% ---- Styling ----
lw = 3;
fs_label = 25;
fs_tick  = 20;
axlw = 1.5;

Fx_on_t1 = 80; Fx_on_t2 = 160;
Fx_bg_color = [1 0 0];
Fx_bg_alpha = 0.1;   % 더 투명하게 원하면 0.15~0.25 추천

%% =========================================================
% FIXED FIGURE GEOMETRY (same for Fig1 and Fig2)
% We'll build 3×3 axes with explicit positions in normalized figure units.
%% =========================================================
figW = 1700; figH = 1150;

% Outer margins (normalized)
L = 0.07; R = 0.03; B = 0.08; T = 0.04;

% Gaps between panels (normalized)  <-- 여기서 "subplot 간격"이 결정됨
gx = 0.035;
gy = 0.045;

% Precompute 3×3 axes positions (row-major, row1=top)
pos = compute_grid_positions(3,3,L,R,B,T,gx,gy);

%% =========================
% FIGURE 1 (3×3): non-TF #1–9
%% =========================
f1 = figure('Color','w','Units','pixels','Position',[50 50 figW figH]);

for k=1:9
    ax = axes('Parent',f1,'Position',pos{k});
    hold(ax,'on'); box(ax,'on'); ax.LineWidth = axlw;

    g = nonTF_idx(k);
    plot(ax, tgrid, X_switch(:,g), 'LineWidth', lw);
    add_time_window_shade(ax, Fx_on_t1, Fx_on_t2, Fx_bg_color, Fx_bg_alpha);
    xline(ax, 80, ':', 'LineWidth', axlw);

    text(ax, 0.02, 0.92, genes(g), 'Units','normalized', ...
        'FontSize', fs_label, 'Interpreter','none', ...
        'HorizontalAlignment','left','VerticalAlignment','top');

    set(ax,'FontSize',fs_tick,'FontName','Times New Roman');
    xticks(ax, xt_all); xlim(ax,[0 160]);
end

exportgraphics(f1, 'Fig9_1_rows1to3.png', 'Resolution', 600);
disp('✅ Saved: Fig9_1_rows1to3.png');

%% =========================
% FIGURE 2 (SAME 3×3 positions): non-TF #10–13 + ΔAUC + heatmap + blanks
%% =========================
f2 = figure('Color','w','Units','pixels','Position',[120 50 figW figH]);

% tiles 1–4: non-TF 10–13
for kk=10:13
    ax = axes('Parent',f2,'Position',pos{kk-9}); % pos{1..4}
    hold(ax,'on'); box(ax,'on'); ax.LineWidth = axlw;

    g = nonTF_idx(kk);
    plot(ax, tgrid, X_switch(:,g), 'LineWidth', lw);
    add_time_window_shade(ax, Fx_on_t1, Fx_on_t2, Fx_bg_color, Fx_bg_alpha);
    xline(ax, 80, ':', 'LineWidth', axlw);

    text(ax, 0.02, 0.92, genes(g), 'Units','normalized', ...
        'FontSize', fs_label, 'Interpreter','none', ...
        'HorizontalAlignment','left','VerticalAlignment','top');

    set(ax,'FontSize',fs_tick,'FontName','Times New Roman');
    xticks(ax, xt_all); xlim(ax,[0 160]);
end

% tile 5: ΔAUC bar  (pos{5})
ax14 = axes('Parent',f2,'Position',pos{7});
hold(ax14,'on'); box(ax14,'on'); ax14.LineWidth = axlw;

bar(ax14, dAUC);
yline(ax14, 0, 'k', 'LineWidth', axlw);
xticks(ax14, 1:N);
xticklabels(ax14, genes_nonTF);
xtickangle(ax14, 45);
ylabel(ax14, '\DeltaAUC_{80-160} (switch - noFx)', 'FontSize', fs_label, 'Interpreter','tex');
set(ax14,'FontSize',fs_tick,'FontName','Times New Roman');

% tile 6: heatmap (pos{6}) + manual colorbar (does NOT resize axes)
ax15 = axes('Parent',f2,'Position',pos{6});
box(ax15,'on'); ax15.LineWidth = axlw;

imagesc(ax15, t_post, 1:N, D.');
set(ax15,'YDir','normal');
hold(ax15,'on');

% Overlay near-zero as gray
tol = 1e-8;
mask0 = abs(D) < tol; mask0 = mask0.';   % (N×T)
zero_gray = 0.88;
rgb = ones(N, numel(t_post), 3) * zero_gray;
h0 = image(ax15, 'XData',[t_post(1) t_post(end)], 'YData',[1 N], 'CData', rgb);
set(h0, 'AlphaData', 0.95 * double(mask0));
uistack(h0,'top');

xticks(ax15,[80 120 160]);
yticks(ax15,1:N);
yticklabels(ax15, genes_nonTF);
xlabel(ax15,'Time (min)','FontSize',fs_label);
set(ax15,'FontSize',fs_tick,'FontName','Times New Roman');

text(ax15, 0.02, 0.98, 'Fx-induced divergence  x^{switch}(t) - x^{noFx}(t)', ...
    'Units','normalized','FontSize',fs_label,'Interpreter','tex', ...
    'HorizontalAlignment','left','VerticalAlignment','top');

% Manual colorbar position (to the right of tile 6, fixed)
cbPos = pos{6};
cbW = 0.015;
cbGap = 0.008;
cbX = cbPos(1) + cbPos(3) + cbGap;
cbY = cbPos(2);
cbH = cbPos(4);

cb = colorbar(ax15,'Position',[cbX cbY cbW cbH]);
cb.FontName = 'Times New Roman';
cb.FontSize = fs_tick;

% --- Gray "0" swatch next to the colorbar (drawn as its own axes) ---
sw_w = cbW * 0.9;
sw_h = cbH * 0.05;
sw_x = cbX + cbW + 0.007;
sw_y = cbY + cbH*0.50 - sw_h/2;

swPos = [sw_x sw_y sw_w sw_h];

% clip to figure bounds so it never exceeds [0,1]
swPos(1) = min(swPos(1), 0.99 - swPos(3));
swPos(2) = max(0.01, min(swPos(2), 0.99 - swPos(4)));

axSw = axes('Parent', f2, 'Units','normalized', 'Position', swPos);
axis(axSw,'off');
rectangle(axSw,'Position',[0 0 1 1], ...
          'FaceColor',[0.5 0.5 0.5], 'EdgeColor','none');

% optional "0" label
text(axSw, 1.2, 0.5, '0', 'Units','normalized', ...
     'VerticalAlignment','middle', ...
     'FontName','Times New Roman', 'FontSize', fs_tick);


% annotation(f2,'rectangle',[sw_x sw_y sw_w sw_h], ...
%     'FaceColor',[zero_gray zero_gray zero_gray], 'EdgeColor','k','LineWidth',1);
% annotation(f2,'textbox',[sw_x+sw_w+0.005 sw_y-0.005 0.04 sw_h*2], ...
%     'String','0','LineStyle','none', ...
%     'FontName','Times New Roman','FontSize',fs_tick, ...
%     'VerticalAlignment','middle','HorizontalAlignment','left');

% tiles 7–9 blanks: just do nothing OR create invisible axes
for b=7:9
    axb = axes('Parent',f2,'Position',pos{b});
    axis(axb,'off');
end

exportgraphics(f2, 'Fig9_2_rows4to5_plus_blank.png', 'Resolution', 600);
disp('✅ Saved: Fig9_2_rows4to5_plus_blank.png');

%% =========================================================
% ---- Local functions ----
%% =========================================================
function pos = compute_grid_positions(nr,nc,L,R,B,T,gx,gy)
% Returns cell array pos{k} (k=1..nr*nc), row-major, row1=top.
W = 1 - L - R;
H = 1 - B - T;
axW = (W - (nc-1)*gx)/nc;
axH = (H - (nr-1)*gy)/nr;

pos = cell(nr*nc,1);
k = 0;
for r = 1:nr
    for c = 1:nc
        k = k + 1;
        x = L + (c-1)*(axW+gx);
        y = B + (nr-r)*(axH+gy);
        pos{k} = [x y axW axH];
    end
end
end

function dx = local_rhs(x, Fxflag, alpha, delta, gamma, W, isTF, K, hill_n)
G = numel(x);
dx = alpha - delta .* x;
if Fxflag==1, dx = dx + gamma; end

tf_idx = find(isTF);
tf_map = zeros(G,1);
tf_map(tf_idx) = 1:numel(tf_idx);

for i = 1:G
    acc = 0;
    for j = 1:G
        if j==i, continue; end
        Wij = W(i,j);
        if Wij==0, continue; end

        if isTF(j)
            kidx = tf_map(j);
            if kidx<1 || kidx>numel(K), continue; end
            Kj = K(kidx);
            xj = max(x(j),0);
            xpow = xj^hill_n;
            den  = (Kj^hill_n) + xpow + 1e-12;
            h    = xpow / den;
            acc  = acc + Wij * h;
        else
            acc  = acc + Wij * x(j);
        end
    end
    dx(i) = dx(i) + acc;
end
end

function add_time_window_shade(ax, x1, x2, faceColor, faceAlpha)
yl = ylim(ax);
p = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], faceColor, ...
    'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
uistack(p,'bottom');
ylim(ax, yl);
end


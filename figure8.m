%% =========================================================
% Supplementary Figure:
% Network-state snapshots (circle layout, Figure 6B style)
% - Edge: fixed inferred network W
% - Activation = green arrow
% - Inhibition = red hammer bar
% - Node: color changes over time
% - TF nodes: binary ON/OFF based on x_j(t) >= K_j
% - If a TF is OFF, outgoing edges from that TF are not drawn
% - Time points: 40, 80, 120, 160 min
% - 0-80 min: IL-1beta only
% - 80-160 min: IL-1beta + Fx
% - Draw one snapshot per figure
%% =========================================================
clear; clc; close all;

%% ---- Global font defaults ----
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultLegendFontName','Times New Roman');
set(groot,'defaultColorbarFontName','Times New Roman');

%% ---- Files ----
fit_file   = 'fit_result.mat';
genes_file = 'genes_used.txt';

assert(exist(fit_file,'file')==2, 'Cannot find %s', fit_file);
assert(exist(genes_file,'file')==2, 'Cannot find %s', genes_file);

S = load(fit_file);

genes = string(readlines(genes_file));
genes = strtrim(genes);
genes = genes(genes ~= "");

%% ---- Required variables from fit_result.mat ----
W      = S.W;
alpha  = S.alpha(:);
delta  = S.delta(:);
gamma  = S.gamma(:);
s      = S.s(:);
isTF   = logical(S.isTF(:));
K      = S.K(:);
hill_n = S.hill_n;

G = size(W,1);
assert(size(W,2)==G, 'W must be square.');
assert(numel(genes)==G, 'genes_used.txt length does not match W size.');

%% ---- Time settings ----
t1 = (0:1:80).';      % IL-1beta only
t2 = (80:1:160).';    % IL-1beta + Fx
tFull = (0:1:160).';

ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',2);

%% ---- Initial condition ----
x0 = max(s, 1e-8);

%% ---- Piecewise simulation ----
rhs_none = @(t,x) local_rhs_snapshot(x, 0, alpha, delta, gamma, W, isTF, K, hill_n);
rhs_fx   = @(t,x) local_rhs_snapshot(x, 1, alpha, delta, gamma, W, isTF, K, hill_n);

[~, X1] = ode15s(rhs_none, t1, x0, ode_opts);
x80 = X1(end,:).';
[~, X2] = ode15s(rhs_fx, t2, x80, ode_opts);

Xfull = [X1; X2(2:end,:)];  % remove duplicated t=80

%% ---- Non-TF normalization: x_i(t) / max_t x_i(t) ----
Xnorm = zeros(size(Xfull));
for i = 1:G
    denom = max(Xfull(:,i));
    if denom < 1e-12
        Xnorm(:,i) = zeros(size(Xfull(:,i)));
    else
        Xnorm(:,i) = Xfull(:,i) ./ denom;
    end
end

%% ---- TF threshold ON/OFF states ----
% K is assumed to correspond to TFs in the order find(isTF)
tf_idx = find(isTF);
nTF = numel(tf_idx);
assert(numel(K) >= nTF, 'K must contain at least one threshold per TF.');

TF_on = false(size(Xfull,1), G);
for kk = 1:nTF
    j = tf_idx(kk);
    TF_on(:,j) = Xfull(:,j) >= K(kk);
end

%% ---- Snapshot times ----
snapTimes  = [40 80 120 160];
snapTitles = { ...
    '40 min (IL-1\beta only)', ...
    '80 min (pre-Fx)', ...
    '120 min (IL-1\beta + Fx)', ...
    '160 min (IL-1\beta + Fx)'};

snapIdx = zeros(size(snapTimes));
for k = 1:numel(snapTimes)
    [~, snapIdx(k)] = min(abs(tFull - snapTimes(k)));
end

%% ---- Fixed network from W ----
thr = 0.01;
W2 = W;
W2(1+G*(0:G-1)) = 0;

mask = ~eye(G) & (abs(W2) > thr);
[ti, sj] = find(mask);         % source j -> target i
weights  = W2(sub2ind([G G], ti, sj));

if isempty(weights)
    error('No edges above threshold %.3g', thr);
end

%% ---- Circle layout ----
theta = linspace(0, 2*pi, G+1).';
theta(end) = [];
R = 2.5;
X = R*cos(theta);
Y = R*sin(theta);
axisScale = 1.2;

%% ---- Visual settings ----
nodeR     = 0.22;
arrowLen  = 0.14;
arrowWid  = 0.09;
barHalf   = 0.09;
gap       = 0.035;
edgeLW    = 1.8;

edgeColorAct = [0.10 0.60 0.20];  % green arrow
edgeColorInh = [0.85 0.15 0.15];  % red hammer bar

fsNode  = 11;
fsTitle = 18;

cmap = parula(256);

%% ---- Draw one snapshot per figure ----
for p = 1:4
    f = figure('Color','w','Position',[100 100 1100 1000]);
    ax = axes('Parent', f);
    hold(ax, 'on');
    axis(ax, 'equal');
    axis(ax, 'off');

    xlim(ax, axisScale * R * [-1 1]);
    ylim(ax, axisScale * R * [-1 1]);

    tIdx = snapIdx(p);

    % ---- Draw nodes first ----
    for i = 1:G
        if isTF(i)
            if TF_on(tIdx,i)
                faceCol = [1.0 0.7 0.7];   % ON: light red
            else
                faceCol = [1 1 1];         % OFF: white
            end
            edgeCol = [0 0 0];
            lw = 3.0;
        else
            cval = Xnorm(tIdx,i);
            faceCol = interp1(linspace(0,1,256), cmap, cval);
            edgeCol = [0.25 0.25 0.30];
            lw = 1.5;
        end

        rectangle(ax, 'Position',[X(i)-nodeR, Y(i)-nodeR, 2*nodeR, 2*nodeR], ...
                     'Curvature',[1 1], ...
                     'FaceColor', faceCol, ...
                     'EdgeColor', edgeCol, ...
                     'LineWidth', lw);

        txt = char(genes(i));
        if numel(txt) > 14
            ksplit = ceil(numel(txt)/2);
            txt = sprintf('%s\n%s', txt(1:ksplit), txt(ksplit+1:end));
        end

        text(ax, X(i), Y(i), txt, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','middle', ...
             'FontWeight','bold', ...
             'FontSize',fsNode, ...
             'Color',[0.10 0.10 0.15], ...
             'Interpreter','none');
    end

    % ---- Draw edges after nodes (on top) ----
    for e = 1:numel(weights)
        i = ti(e);
        j = sj(e);
        wij = weights(e);

        if isTF(j) && ~TF_on(tIdx,j)
            continue;
        end

        P = [X(j), Y(j)];
        Q = [X(i), Y(i)];
        v = Q - P;
        L = hypot(v(1), v(2));
        if L < eps
            continue;
        end

        u = v / L;
        n = [-u(2), u(1)];

        startPt = P + nodeR * u;

        if wij > 0
            % green arrow
            tip = Q - (nodeR + gap) * u;
            base = tip - arrowLen * u;

            line(ax, [startPt(1), base(1)], [startPt(2), base(2)], ...
                 'Color', edgeColorAct, 'LineWidth', edgeLW);

            Lft = base + 0.5 * arrowWid * n;
            Rgt = base - 0.5 * arrowWid * n;
            patch(ax, [tip(1), Lft(1), Rgt(1)], [tip(2), Lft(2), Rgt(2)], ...
                  edgeColorAct, 'EdgeColor', 'none');

        else
            % red hammer bar
            barCenter = Q - (nodeR + gap) * u;
            stemEnd   = barCenter - 0.02 * u;

            line(ax, [startPt(1), stemEnd(1)], [startPt(2), stemEnd(2)], ...
                 'Color', edgeColorInh, 'LineWidth', edgeLW);

            B1 = barCenter + barHalf * n;
            B2 = barCenter - barHalf * n;
            line(ax, [B1(1), B2(1)], [B1(2), B2(2)], ...
                 'Color', edgeColorInh, 'LineWidth', edgeLW + 0.8);
        end
    end

    % ---- Title ----
    title(ax, snapTitles{p}, 'FontWeight','bold', 'FontSize',fsTitle);

    % ---- Colorbar (for non-TF scale) ----
    colormap(ax, parula);
    cb = colorbar(ax);
    cb.FontName = 'Times New Roman';
    cb.FontSize = 13;
    cb.Label.String = {'Node color'; 'non-TF: normalized activity'; 'TF: OFF (white) / ON (light red)'};
    cb.Label.FontSize = 13;

    % ---- Save ----
    outname = sprintf('NetworkSnapshot_%dmin.png', snapTimes(p));
    exportgraphics(f, outname, 'Resolution', 600);
    disp(['Saved: ' outname]);

    % Optional: save pdf too
    % exportgraphics(f, sprintf('NetworkSnapshot_%dmin.pdf', snapTimes(p)), 'ContentType','vector');

    % close(f);
end

%% =========================================================
% Local RHS
%% =========================================================
function dx = local_rhs_snapshot(x, Fxflag, alpha, delta, gamma, W, isTF, K, hill_n)
G = numel(x);
dx = alpha - delta .* x;

if Fxflag == 1
    dx = dx + gamma;
end

tf_idx = find(isTF);
tf_map = zeros(G,1);
tf_map(tf_idx) = 1:numel(tf_idx);

for i = 1:G
    acc = 0;
    for j = 1:G
        if i == j
            continue;
        end

        Wij = W(i,j);
        if Wij == 0
            continue;
        end

        if isTF(j)
            kidx = tf_map(j);
            Kj = K(kidx);

            xj = max(x(j), 0);
            xpow = xj^hill_n;
            h = xpow / (Kj^hill_n + xpow + 1e-12);

            acc = acc + Wij * h;
        else
            acc = acc + Wij * x(j);
        end
    end
    dx(i) = dx(i) + acc;
end
end
function plot_double_lesion_geneheatmap(cfg)
%PLOT_DOUBLE_LESION_GENEHEATMAP
% Supplementary full gene-wise double-lesion interaction heatmap

%% ---------------- defaults ----------------
assert(isfield(cfg,'outDir') && ischar(cfg.outDir), 'cfg.outDir must be set.');

if ~isfield(cfg,'savePrefix') || isempty(cfg.savePrefix)
    cfg.savePrefix = 'Fig_double_lesion_geneheatmap_full';
end
if ~isfield(cfg,'useAsinh') || isempty(cfg.useAsinh)
    cfg.useAsinh = true;
end
if ~isfield(cfg,'kFrac') || isempty(cfg.kFrac)
    cfg.kFrac = 0.05;
end
if ~isfield(cfg,'fontName') || isempty(cfg.fontName)
    cfg.fontName = 'Times New Roman';
end
if ~isfield(cfg,'fsLabel') || isempty(cfg.fsLabel)
    cfg.fsLabel = 13;
end
if ~isfield(cfg,'fsTick') || isempty(cfg.fsTick)
    cfg.fsTick = 10;
end
if ~isfield(cfg,'sortPairs') || isempty(cfg.sortPairs)
    cfg.sortPairs = true;
end
if ~isfield(cfg,'sortGenes') || isempty(cfg.sortGenes)
    cfg.sortGenes = true;
end

%% ---------------- load ----------------
S = load(fullfile(cfg.outDir,'double_lesion_results.mat'));
assert(isfield(S,'results2'), 'double_lesion_results.mat must contain results2');
R2 = S.results2;

req = {'I_gene_pair','idxOut','genes','pairs','lesionTargets'};
for k = 1:numel(req)
    assert(isfield(R2,req{k}), 'results2 missing field: %s', req{k});
end

I = R2.I_gene_pair;
idxOut = R2.idxOut(:);
genes = R2.genes;
pairs = R2.pairs;
lesionTargets = R2.lesionTargets(:);

[nOut,nPairs] = size(I);

%% ---------------- labels ----------------
pairLabels = cell(nPairs,1);
for p = 1:nPairs
    a = pairs(p,1);
    b = pairs(p,2);
    j = lesionTargets(a);
    k = lesionTargets(b);
    pairLabels{p} = sprintf('%s × %s', genes{j}, genes{k});
end
geneLabels = genes(idxOut);

%% ---------------- sorting ----------------
if cfg.sortPairs
    if isfield(R2,'I_meanAbs') && numel(R2.I_meanAbs)==nPairs
        pairScore = R2.I_meanAbs(:);
    else
        pairScore = mean(abs(I),1).';
    end
    [~, ordP] = sort(pairScore, 'descend');
    I = I(:, ordP);
    pairLabels = pairLabels(ordP);
end

if cfg.sortGenes
    geneScore = mean(abs(I),2);
    [~, ordG] = sort(geneScore, 'descend');
    I = I(ordG,:);
    geneLabels = geneLabels(ordG);
end

%% ---------------- scaling ----------------
absI = abs(I(:));
absI = absI(isfinite(absI));
if isempty(absI)
    clim = 1;
else
    clim = prctile(absI, 98);
    if clim == 0
        clim = max(absI);
        if clim == 0, clim = 1; end
    end
end

if cfg.useAsinh
    k = cfg.kFrac * clim;
    Iviz = asinh(I / k);
    climT = asinh(clim / k);
    cLimUse = [-climT, climT];
    cbTickPos = [-climT, 0, climT];
    cbTickLab = {sprintf('-%0.2g',clim), '0', sprintf('%0.2g',clim)};
else
    Iviz = I;
    cLimUse = [-clim, clim];
    cbTickPos = [-clim, 0, clim];
    cbTickLab = {sprintf('-%0.2g',clim), '0', sprintf('%0.2g',clim)};
end

%% ---------------- plot ----------------
f = figure('Color','w','Position',[60 60 1550 780]);
ax = axes('Parent',f);

imagesc(ax, Iviz);
set(ax,'YDir','normal');
box(ax,'on');
set(ax,'FontName',cfg.fontName,'FontSize',cfg.fsTick);
caxis(ax, cLimUse);
colormap(ax, make_diverging_cmap(256, 0.97));

xlabel(ax, 'TF pair (simultaneous lesion)', 'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
ylabel(ax, 'Downstream gene', 'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
title(ax, 'Gene-wise double-lesion interaction  I_{i,kl}  (full)', ...
    'FontName', cfg.fontName, 'FontSize', cfg.fsLabel+1, 'FontWeight','bold');

xticks(ax, 1:nPairs);
xticklabels(ax, pairLabels);
xtickangle(ax, 45);

yticks(ax, 1:numel(geneLabels));
yticklabels(ax, geneLabels);

ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridColor = [0.88 0.88 0.88];
ax.GridAlpha = 0.35;
ax.TickLength = [0 0];

cb = colorbar(ax);
cb.FontName = cfg.fontName;
cb.FontSize = cfg.fsTick;
cb.Ticks = cbTickPos;
cb.TickLabels = cbTickLab;
cb.Label.String = 'Interaction  I_{i,kl}  (synergy + / compensation -)';
cb.Label.FontSize = cfg.fsLabel;

outPDF = fullfile(cfg.outDir, [cfg.savePrefix '.pdf']);
outPNG = fullfile(cfg.outDir, [cfg.savePrefix '.png']);
exportgraphics(f, outPDF, 'ContentType','vector');
exportgraphics(f, outPNG, 'Resolution',300);

fprintf('Saved:\n  %s\n  %s\n', outPDF, outPNG);
end

%% =========================================================
function cmap = make_diverging_cmap(n, midGray)
n2 = floor(n/2);
c1 = [linspace(0, midGray, n2)' ...
      linspace(0, midGray, n2)' ...
      linspace(1, midGray, n2)'];
c2 = [linspace(midGray, 1, n-n2)' ...
      linspace(midGray, 0, n-n2)' ...
      linspace(midGray, 0, n-n2)'];
cmap = [c1; c2];
end
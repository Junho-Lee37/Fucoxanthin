function plot_lesion_multipanel(cfg)
%PLOT_LESION_MULTIPANEL
% Figure 9 in 2x2 format:
%   A: TF causal contribution
%   B: single-lesion heatmap
%   C: direct gamma vs network sensitivity
%   D: double-lesion interaction heatmap
%
% mode:
%   cfg.mode = 'top'  -> D panel uses top TF pairs only
%   cfg.mode = 'full' -> D panel uses full TF pairs

%% ---------------- defaults ----------------
assert(isfield(cfg,'outDir') && ischar(cfg.outDir), 'cfg.outDir must be set.');
assert(isfield(cfg,'mode') && (strcmp(cfg.mode,'top') || strcmp(cfg.mode,'full')), ...
    'cfg.mode must be ''top'' or ''full''');

if ~isfield(cfg,'savePrefix') || isempty(cfg.savePrefix)
    cfg.savePrefix = 'Fig_lesion';
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
if ~isfield(cfg,'kFrac') || isempty(cfg.kFrac)
    cfg.kFrac = 0.05;
end
if ~isfield(cfg,'topNPairsMain') || isempty(cfg.topNPairsMain)
    cfg.topNPairsMain = 8;
end
if ~isfield(cfg,'sortPairs') || isempty(cfg.sortPairs)
    cfg.sortPairs = true;
end
if ~isfield(cfg,'sortGenes') || isempty(cfg.sortGenes)
    cfg.sortGenes = true;
end

%% ---------------- load data ----------------
B = load(fullfile(cfg.outDir,'baseline.mat'));
L = load(fullfile(cfg.outDir,'lesion_results.mat'));
D = load(fullfile(cfg.outDir,'double_lesion_results.mat'));

assert(isfield(B,'deltaAUC') && isfield(B,'genes') && isfield(B,'idxOut'), ...
    'baseline.mat missing required fields');
assert(isfield(L,'results'), 'lesion_results.mat missing results');
assert(isfield(D,'results2'), 'double_lesion_results.mat missing results2');

R  = L.results;
R2 = D.results2;

genes  = B.genes;
idxOut = B.idxOut(:);
dBase  = B.deltaAUC(:);

%% =========================================================
% Panel A/B/C data: single lesion
%% =========================================================
lesions = R.lesionTargets(:);
nLes    = numel(lesions);
nOut    = numel(idxOut);

dLes = zeros(nOut, nLes);
for m = 1:nLes
    dLes(:,m) = R.perLesion(m).deltaAUC(:);
end

% ΔΔAUC = lesion response - wild-type response
dDiff = dLes - dBase;                  % nOut x nLes
scoreTF = mean(abs(dDiff),1);          % 1 x nLes

% sort TFs by mean absolute lesion effect
[scoreS, ordTF] = sort(scoreTF(:), 'descend');
tfIdxSorted = lesions(ordTF);
tfNamesS = genes(tfIdxSorted);
Hsingle  = dDiff(:, ordTF);            % gene x TF, sorted

% gamma for panel C
if isfield(B,'gamma')
    gam = B.gamma(:);
elseif isfield(cfg,'fitMatFile') && ~isempty(cfg.fitMatFile)
    Sfit = load(cfg.fitMatFile);
    assert(isfield(Sfit,'gamma'), 'fit_result.mat missing gamma');
    gam = Sfit.gamma(:);
else
    error('gamma not found in baseline.mat or cfg.fitMatFile');
end
geneNet = mean(abs(dDiff),2);

%% =========================================================
% Panel D data: double lesion
%% =========================================================
Iall = R2.I_gene_pair;              % nOut x nPairs
pairsAll = R2.pairs;
lesionTargets = R2.lesionTargets(:);

[nOut2, nPairsAll] = size(Iall);
assert(nOut2 == nOut, 'Mismatch between single/double-lesion downstream genes.');

pairLabelsAll = cell(nPairsAll,1);
for p = 1:nPairsAll
    a = pairsAll(p,1);
    b = pairsAll(p,2);
    j = lesionTargets(a);
    k = lesionTargets(b);
    pairLabelsAll{p} = sprintf('%s × %s', genes{j}, genes{k});
end

% sort pairs by interaction strength
if cfg.sortPairs
    if isfield(R2,'I_meanAbs') && numel(R2.I_meanAbs)==nPairsAll
        pairScoreAll = R2.I_meanAbs(:);
    else
        pairScoreAll = mean(abs(Iall),1).';
    end
    [~, ordP] = sort(pairScoreAll, 'descend');
    I_D = Iall(:, ordP);
    pairLabels_D = pairLabelsAll(ordP);
else
    I_D = Iall;
    pairLabels_D = pairLabelsAll;
end

% mode-specific selection for panel D
if strcmp(cfg.mode,'top')
    nKeep = min(cfg.topNPairsMain, size(I_D,2));
    I_D = I_D(:, 1:nKeep);
    pairLabels_D = pairLabels_D(1:nKeep);
    titleD = sprintf('(D) Double-lesion interaction (top %d pairs)', nKeep);
else
    titleD = '(D) Double-lesion interaction (full)';
end

% sort downstream genes for panel D only
geneLabels_D = genes(idxOut);
if cfg.sortGenes
    geneScoreD = mean(abs(I_D),2);
    [~, ordG_D] = sort(geneScoreD, 'descend');
    I_D = I_D(ordG_D,:);
    geneLabels_D = geneLabels_D(ordG_D);
end

% scale panel D
absI = abs(I_D(:));
absI = absI(isfinite(absI));
if isempty(absI)
    climD = 1;
else
    climD = prctile(absI, 98);
    if climD == 0
        climD = max(absI);
        if climD == 0, climD = 1; end
    end
end

kD = cfg.kFrac * climD;
Iviz_D = asinh(I_D / kD);
climTD = asinh(climD / kD);

%% =========================================================
% Figure layout: 2 x 2
%% =========================================================
f = figure('Color','w','Position',[70 70 1450 1050]);

%% ---------------- A ----------------
axA = subplot(2,2,1);
bar(axA, scoreS, 'FaceColor',[0.20 0.40 0.70], 'EdgeColor','none');
box(axA,'on'); grid(axA,'on');
set(axA,'FontName',cfg.fontName,'FontSize',cfg.fsTick);

xticks(axA, 1:nLes);
xticklabels(axA, tfNamesS);
xtickangle(axA, 45);

ymax = max(scoreS);
if ymax <= 0, ymax = 1; end
ylim(axA, [0, 1.10*ymax]);
xlim(axA, [0.5, nLes+0.5]);

ylabel(axA, 'Mean |\Delta\DeltaAUC| across downstream genes', ...
    'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
title(axA, '(A) TF causal contribution', ...
    'FontName', cfg.fontName, 'FontSize', cfg.fsLabel+2, 'FontWeight','bold');

%% ---------------- B ----------------
axB = subplot(2,2,2);

absH = abs(Hsingle(:));
absH = absH(isfinite(absH));
if isempty(absH)
    climB = 1;
else
    climB = prctile(absH, 98);
    if climB == 0
        climB = max(absH);
        if climB == 0, climB = 1; end
    end
end

kB = cfg.kFrac * climB;
Hviz = asinh(Hsingle / kB);
climTB = asinh(climB / kB);

imagesc(axB, Hviz);
set(axB,'YDir','normal');
box(axB,'on');
set(axB,'FontName',cfg.fontName,'FontSize',cfg.fsTick);
caxis(axB, [-climTB, climTB]);
colormap(axB, make_diverging_cmap(256, 0.95));

xlabel(axB, 'Lesioned TF (ranked)', 'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
ylabel(axB, 'Downstream gene', 'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
title(axB, '(B) \Delta\DeltaAUC by single TF lesion', ...
    'FontName', cfg.fontName, 'FontSize', cfg.fsLabel+2, 'FontWeight','bold');

xticks(axB, 1:nLes);
xticklabels(axB, tfNamesS);
xtickangle(axB, 45);

if nOut <= 25
    yticks(axB, 1:nOut);
    yticklabels(axB, genes(idxOut));
else
    yt = unique(round(linspace(1,nOut,12)));
    yticks(axB, yt);
    yticklabels(axB, genes(idxOut(yt)));
end

cbB = colorbar(axB);
cbB.FontName = cfg.fontName;
cbB.FontSize = cfg.fsTick;
cbB.Ticks = [-climTB, 0, climTB];
cbB.TickLabels = {sprintf('-%0.2g',climB), '0', sprintf('%0.2g',climB)};
cbB.Label.String = '\Delta\DeltaAUC';
cbB.Label.FontSize = cfg.fsLabel-1;

%% ---------------- C ----------------
axC = subplot(2,2,3);
box(axC,'on'); grid(axC,'on'); hold(axC,'on');
set(axC,'FontName',cfg.fontName,'FontSize',cfg.fsTick);

scatter(axC, gam(idxOut), geneNet, 60, 'filled', ...
    'MarkerFaceAlpha',0.78, 'MarkerEdgeColor','k');

xlabel(axC, 'Direct Fx effect  \gamma_i', 'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
ylabel(axC, 'Mean network sensitivity', 'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
title(axC, '(C) Direct vs network-mediated contribution', ...
    'FontName', cfg.fontName, 'FontSize', cfg.fsLabel+2, 'FontWeight','bold');

%% ---------------- D ----------------
axD = subplot(2,2,4);

imagesc(axD, Iviz_D);
set(axD,'YDir','normal');
box(axD,'on');
set(axD,'FontName',cfg.fontName,'FontSize',cfg.fsTick);
caxis(axD, [-climTD, climTD]);
colormap(axD, make_diverging_cmap(256, 0.97));

xlabel(axD, 'TF pair (simultaneous lesion)', 'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
ylabel(axD, 'Downstream gene', 'FontName', cfg.fontName, 'FontSize', cfg.fsLabel);
title(axD, titleD, ...
    'FontName', cfg.fontName, 'FontSize', cfg.fsLabel+2, 'FontWeight','bold');

xticks(axD, 1:numel(pairLabels_D));
xticklabels(axD, pairLabels_D);
xtickangle(axD, 45);

yticks(axD, 1:numel(geneLabels_D));
yticklabels(axD, geneLabels_D);

axD.XGrid = 'on';
axD.YGrid = 'on';
axD.GridColor = [0.88 0.88 0.88];
axD.GridAlpha = 0.35;
axD.TickLength = [0 0];

cbD = colorbar(axD);
cbD.FontName = cfg.fontName;
cbD.FontSize = cfg.fsTick;
cbD.Ticks = [-climTD, 0, climTD];
cbD.TickLabels = {sprintf('-%0.2g',climD), '0', sprintf('%0.2g',climD)};
cbD.Label.String = 'Interaction  I_{i,kl}';
cbD.Label.FontSize = cfg.fsLabel-1;

%% ---------------- save ----------------
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
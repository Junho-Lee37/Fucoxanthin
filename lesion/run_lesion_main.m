clc; clear;

cfg.fitMatFile = 'fit_result.mat';
cfg.genesFile  = 'genes_used.txt';
cfg.outDir     = 'lesion_results';

cfg.tspan      = [0 160];
cfg.aucWindow  = [80 160];
cfg.nGrid      = 1601;
cfg.fxOn       = 80;

cfg.lesionTargets = [];
cfg.useNonTFOnly  = true;
cfg.doBaseline    = true;

cfg.fontName      = 'Times New Roman';
cfg.fsLabel       = 13;
cfg.fsTick        = 10;
cfg.kFrac         = 0.01; % grey

cfg.topNPairsMain = 8;
cfg.sortPairs     = true;
cfg.sortGenes     = true;

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

run_network_lesion(cfg);
run_double_lesion(cfg);

%% -------- Figure 1: main (TOP D) --------
cfg.mode = 'top';
cfg.savePrefix = 'Fig_lesion_main';
plot_lesion_multipanel(cfg);

%% -------- Figure 2: supplementary (FULL D) --------
cfg.mode = 'full';
cfg.savePrefix = 'Fig_lesion_full';
plot_lesion_multipanel(cfg);

fprintf('\nDone. Generated:\n');
fprintf('  Fig_lesion_main (top D)\n');
fprintf('  Fig_lesion_full (full D)\n');
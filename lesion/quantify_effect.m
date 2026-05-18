function q = quantify_effect(t_none, x_none, t_fx, x_fx, aucWindow, idxOut, nGrid)
%QUANTIFY_EFFECT
% deltaAUC = AUC(Fx) - AUC(None) for each gene (restricted to idxOut)

if nargin < 7 || isempty(nGrid)
    nGrid = 401;
end

tGrid = linspace(aucWindow(1), aucWindow(2), nGrid);

Xn = interp1(t_none, x_none, tGrid, 'pchip');
Xf = interp1(t_fx,   x_fx,   tGrid, 'pchip');

AUC_none = trapz(tGrid, Xn);
AUC_fx   = trapz(tGrid, Xf);

deltaAUC = AUC_fx - AUC_none;

q.AUC_none = AUC_none(idxOut);
q.AUC_fx   = AUC_fx(idxOut);
q.deltaAUC = deltaAUC(idxOut);
end

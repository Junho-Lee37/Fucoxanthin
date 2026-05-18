function dx = grn_rhs(~, x, theta, Fxflag)
%GRN_RHS
% dx = alpha - delta.*x + W*f(x) + gamma*Fxflag

G = numel(x);
dx = zeros(G,1);

% Expand K to length G if needed
tfIdx = find(theta.isTF);
if numel(theta.K) == numel(tfIdx)
    Kfull = nan(G,1);
    Kfull(tfIdx) = theta.K(:);
elseif numel(theta.K) == G
    Kfull = theta.K(:);
else
    error('K length mismatch: len(K)=%d, nTF=%d, G=%d', numel(theta.K), numel(tfIdx), G);
end

n = theta.hill_n;

% f_j(x_j)
fj = zeros(G,1);
for j = 1:G
    if theta.isTF(j)
        Kj = Kfull(j);
        if ~isfinite(Kj) || Kj <= 0
            fj(j) = x(j); % safe fallback
        else
            xn = x(j)^n;
            fj(j) = xn / (xn + Kj^n);
        end
    else
        fj(j) = x(j);
    end
end

coupling = theta.W * fj;

dx = theta.alpha - theta.delta .* x + coupling + theta.gamma .* Fxflag;
end

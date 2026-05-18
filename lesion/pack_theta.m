function theta = pack_theta(alpha, delta, gamma, W, s, K, isTF, hill_n)
%PACK_THETA pack parameters into a struct for simulation
theta.alpha  = alpha(:);
theta.delta  = delta(:);
theta.gamma  = gamma(:);
theta.W      = W;
theta.s      = s(:);
theta.K      = K(:);
theta.isTF   = logical(isTF(:));
theta.hill_n = hill_n;
end

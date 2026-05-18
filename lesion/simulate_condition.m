function [t, x] = simulate_condition(theta, condName, tspan)
%SIMULATE_CONDITION
% condName: 'None' or 'Fx'

switch condName
    case 'None'
        Fxflag = 0;
    case 'Fx'
        Fxflag = 1;
    otherwise
        error('Unknown condition: %s', condName);
end

% Initial condition: your Stage 3 is x(0)=s
x0 = theta.s;

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1.0);

[t, x] = ode15s(@(t,x) grn_rhs(t, x, theta, Fxflag), tspan, x0, opts);
end

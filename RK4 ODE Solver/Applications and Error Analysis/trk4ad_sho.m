format long

function dydt = harmonic(t, y)
% Stores harmonic oscillator ODE: x'' = -x 
% Inputs
% y: 2D column vector in the form [x; x']
% t: Current time, unused but needed for rk4step function signature
%
% Output
% dydt: Derivative of y, 2D column vector in form [x', x'']

    dydt = [y(2); -y(1)];

end

% Tolerances
reltols = [1.0e-5, 1.0e-7, 1.0e-9, 1.0e-11];

% Preset timespan for all iterations
tspan = linspace(0, 3*pi, 65);

% Exact solution to the H-O with given initial conds
exact_sol = sin(tspan);

% Initial conds
y0 = [0; 1];

% Calculate and plot for each tolerance
hold on;
for reltol = reltols
    % Compute solution
    [tout, yout] = rk4ad(@harmonic, tspan, reltol, y0);

    % Calculate error
    err = exact_sol - yout(1, :);

    % Plot error
    plot(tout, err, "DisplayName", string(reltol)+" Tolerance")
end
legend("Location","southwest");
xlabel("Time");
ylabel("Errors");
title("Errors vs Time For H-O (rk4ad)");
drawnow;

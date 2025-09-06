format long

% Using simple harmonic oscillator as test ODE system, will be starting at
% x(t0) = 1, x'(t0) = v(t0) = 0, so answer will be cosine function 

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

dt_values = [0.5, 0.5^2, 0.5^3, 0.5^4, 0.5^5, 0.5^6, 0.5^7, 0.5^8];
errors = zeros(size(dt_values));

for i = 1:length(dt_values)

    dt = dt_values(i);

    % x = y0(1), x' = v = y0(2)
    y0 = [1; 0];
    
    t0 = 0;
    
    y_step = rk4step(@harmonic, t0, dt, y0);
    
    % Exact solution with given boundary conds are x(t) = cos(t), 
    % x'(t) = -sin(t)
    y_exact = [cos(t0 + dt); -sin(t0 + dt)];
    
    % error should be on order dt^5
    error = norm(y_exact - y_step);
    
    errors(i) = error;
end

ratios = errors(1:end-1) ./ errors(2:end);
disp(ratios);
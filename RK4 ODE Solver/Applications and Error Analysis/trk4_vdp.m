format long

global a;
a = 5;

function dydt = vanderpol(t, y)
% Stores Van der Pol oscillator ODE: x'' + a(x^2 - 1)x' + x = 0 
% Also written as x'' = -a(x^2 - 1)x' - x
% Inputs
% y: 2D column vector in the form [x; x']
% t: Current time, unused but needed for rk4step function signature
%
% Output
% dydt: Derivative of y, 2D column vector in form [x', x'']
    global a;
    dydt = [y(2); -a*(y(1)^2 - 1) * y(2) - y(1)];

end

% Setup time domain
level = 12;
tmax = 100;
nt = 2^(level) + 1;
dt = tmax / 2^(level);
tspan = linspace(0, tmax, nt);

% Initial conds y0 = [x, x']
y0 = [1; -6];

% Run RK4 algorithm
[tout, yout] = rk4(@vanderpol, tspan, y0);

% Plot x(t) vs t
figure;
plot(tout, yout(1, :));
xlabel("$t$",'Interpreter','latex');
ylabel("$x(t)$",'Interpreter','latex');
title("Van der Pol Oscillator Position vs Time (rk4)");
grid on;

% Plot x'(t) vs x(t)
figure;
plot(yout(1, :), yout(2, :));
xlabel("$x(t)$",'Interpreter','latex');
ylabel("$\frac{dx}{dt}(t)$",'Interpreter','latex');
title("VDP Oscillator Phase-Space Velocity vs Position (rk4)",'Interpreter','latex');
grid on;

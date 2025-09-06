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

tmax = 3*pi;

y0 = [0; 1];

% Number points in time for levels 6, 7, 8
nt6 = 2^(6) + 1;
nt7 = 2^(7) + 1;
nt8 = 2^(8) + 1;

% delta time for levels 6, 7, 8
dt6 = tmax / 2^(6);
dt6 = tmax / 2^(7);
dt6 = tmax / 2^(8);

% Time domain values for levels 6, 7, 8
tspan6 = linspace(0, tmax, nt6);
tspan7 = linspace(0, tmax, nt7);
tspan8 = linspace(0, tmax, nt8);

% Solving harmonic oscillator for levels 6, 7, 8
[tout6, yout6] = rk4(@harmonic, tspan6, y0);
[tout7, yout7] = rk4(@harmonic, tspan7, y0);
[tout8, yout8] = rk4(@harmonic, tspan8, y0);

% Finding errors for levels 6, 7, 8
err6 = sin(tspan6) - yout6(1, :);
err7 = sin(tspan7) - yout7(1, :);
err8 = sin(tspan8) - yout8(1, :);


% Plot the scaled errors
rho = 2;
figure;
hold on;
plot(tout6, err6, 'DisplayName', 'Error (Level 6)');
plot(tout7, rho^4 * err7, 'DisplayName', '\rho^4 Error (Level 7)');
plot(tout8, rho^8 * err8, 'DisplayName', '\rho^8 Error (Level 8)');
xlabel("Time");
ylabel("Scaled Errors");
title("Scaled Errors vs Time");
legend("Location", "southwest");
drawnow;
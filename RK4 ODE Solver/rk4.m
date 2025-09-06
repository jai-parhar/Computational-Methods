function [tout, yout] = rk4(fcn, tspan, y0)
% Performs 4th order Runge-Kutta algorithm to solve system of ODEs
% Inputs
% fcn: Function handle for right hand sides of ODEs (returns
% length-n column vector)
% tspan: Vector of output times (length nout).
% y0: Initial values (length-n column vector).
%
% Outputs
% tout: Vector of output times.
% yout: Output values (nout x n array. The ith column of yout
% contains the nout values of the ith dependent variable).

% Find nout and n from inputted vectors
nout = size(tspan, 2);
n = size(y0, 1);

% Initialize output values to be all 0
yout = zeros(n, nout);

% Start at first timeloop, go to second last, as this calculates i+1th
% timestep
for i = 1:nout-1
    % Calculate timestep from input timespans
    dt = tspan(i+1) - tspan(i);
    
    % Store ith timestep
    yout(:, i) = y0;

    % Calculate next step and set to current step
    y_next = rk4step(fcn, tspan(i), dt, y0);
    y0 = y_next;

end

% Store last timestep, not done in loop
yout(:, end) = y0;

% Copies tspan to output times
tout = tspan;

end

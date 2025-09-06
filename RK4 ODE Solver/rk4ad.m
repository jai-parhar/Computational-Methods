function [tout yout] = rk4ad(fcn, tspan, reltol, y0)
% Performs 4th order Runge-Kutta algorithm to solve system of ODEs with
% variable timesteps
% Inputs
% fcn: Function handle for right hand sides of ODEs (returns
% length-n column vector)
% tspan: Vector of output times (length nout vector).
% reltol: Relative tolerance parameter.
% y0: Initial values (length-n column vector).
%
% Outputs
% tout: Output times (length-nout column vector, elements
% identical to tspan).
% yout: Output values (nout x n array. The ith column of yout
% contains the nout values of the ith dependent variable)

% Find nout and n from inputted vectors
nout = size(tspan, 2);
n = size(y0, 1);

% Initialize output values to be all 0
yout = zeros(n, nout);

% Start at first timeloop, go to second last, as this calculates i+1th
% timestep
for i = 1:nout-1
    % Calculate coarse and fine timesteps from input timespans
    dt = tspan(i+1) - tspan(i);
    
    % Store ith timestep
    yout(:, i) = y0;

    % Calculate next coarse and fine step and set to current step
    y_next_C = rk4step(fcn, tspan(i), dt, y0);
    y_next_F_middle = rk4step(fcn, tspan(i), dt/2, y0);
    y_next_F = rk4step(fcn, tspan(i) + dt/2, dt/2, y_next_F_middle);

    % Calculate error between coarse and fine steps
    err_C = (16/15)*norm(y_next_C - y_next_F);

    % Idea to find timestep to use: Take error computed so far, find k(t0)
    % from this, and divide into substeps based on number of times new dt
    % will need to fit in between, then set y0 to value found at end of
    % substep iterations
    
    err_rel = err_C/norm(y_next_F);

    dt_min = 1e-4;

    % Check to see if we need to perform substeps
    if reltol < err_rel
        
        % Calculate the ideal timestep to use that puts us at the threshold
        % Done by scaling dt, then taking max between that and min val
        dt_ideal = max(0.9 * dt * nthroot(reltol/err_rel, 5), dt_min);

        if dt_ideal == 1e-4
            warning("Minimum stepsize of 1e-4 reached or exceeded, using 1e-4 for iteration");
        end

        % Calculate the number of substeps we must take at this new
        % timestep, rounding up
        n_substeps = ceil(dt/dt_ideal) + 1;
        
        % Find times for each substep
        t_substeps = linspace(tspan(i), tspan(i+1), n_substeps);

        % Find what dt must be for each substep
        dt_substep = t_substeps(2) - t_substeps(1);

        dt = dt_substep;

        % SUBSTEP FLOOR
        if dt_substep < dt_min
            n_substeps = floor(dt/dt_min);
            t_substeps = linspace(tspan(i), tspan(i+1), n_substeps);
            dt_substep = t_substeps(2) - t_substeps(1);
        end

        % Initialize substep y
        y_curr_substep = y0;
        for t = t_substeps(1:end-1) 
        % Skip last to make sure not doubling up, dont need to perform
        % iteration at last time anyways because that puts us into 
        % last time + dt_new
            % Perform substep iteration
            %t
            y_next_substep = rk4step(fcn, t, dt_substep, y_curr_substep);
            y_curr_substep = y_next_substep;
        end

        % Store final result
        y0 = y_curr_substep;

    else
        % What we had was good enough, store it
        y0 = y_next_C;
    end

end

% Store last timestep, not done in loop
yout(:, end) = y0;

% Copies tspan to output times
tout = tspan;

end
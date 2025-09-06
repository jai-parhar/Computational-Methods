
function x = hybrid(f, dfdx, xmin, xmax, tol1, tol2)
% f: Function whose root is sought (takes one argument).
% dfdx: Derivative function (takes one argument).
% xmin: Initial bracket minimum.
% xmax: Initial bracket maximum.
% tol1: Relative convergence criterion for bisection.
% tol2: Relative convergence criterion for Newton iteration.
% Output => x: Estimate of root.
    
    fmin = f(xmin);
    fmax = f(xmax);

    bi_max_iterations = 50;

    % Check to make sure we have a root in this interval
    if ~(fmin*fmax < 0)
        error("It is not guaranteed that a root exists in this interval")
    end
    
    bi_converged = false;
    bi_iteration = 0;
    
    % Bisection
    while ~bi_converged

        xmid = (xmax + xmin)/2;
        fmid = f(xmid);

        if fmid == 0.0   % Root found
            bi_converged = true;
        elseif fmin * fmid < 0   % Root on left side of xmid
            xmax = xmid;
        else   % Only option left is that root on right side of xmid
            xmin = xmid;
            fmin = fmid;   % Set new value of fmin for the next iteration
        end
        
        if abs(xmax - xmin)/abs(xmid) < tol1   % Convergence test
            bi_converged = true;
        end

        if bi_iteration > bi_max_iterations   % Check to see if too many iterations have occured
            warning("Exceeded max bisection iterations, breaking with current values")
            break
        end

        bi_iteration = bi_iteration + 1;
    end

    root = xmid;   % We start with our current root bisection for Newton's method
    
    ne_max_iterations = 500;
    ne_iteration = 0;
    ne_converged = false;
    
    curr_x = root;
    while ~ne_converged
        
        next_x = curr_x - (f(curr_x)/dfdx(curr_x));

        if abs(next_x - curr_x)/abs(next_x) < tol2   % Convergence test
            ne_converged = true;
        end

        curr_x = next_x;

        if ne_iteration > ne_max_iterations   % Check to see if too many iterations have occured
            warning("Exceeded max Newton's method iterations, breaking with current values")
            break
        end
        
        ne_iteration = ne_iteration + 1;

    end
    
    x = curr_x;
end


function f_x = f(x)
    f_x = 512*x.^10 - 5120*x.^9 + 21760*x.^8 - 51200*x.^7 + 72800*x.^6 - 64064*x.^5 + ...
        34320*x.^4 - 10560*x.^3 + 1650*x.^2 - 100*x + 1;
end

function df_dx = dfdx(x)
    df_dx = 5120*x.^9 - 46080*x.^8 + 174080*x.^7 - 358400*x.^6 + 436800*x.^5 - 320320*x.^4 + ...
        137280*x.^3 - 31680*x.^2 + 3300*x - 100;
end


x_min = 0;
x_max = 2;
x_vals = linspace(x_min, x_max, 5000);

% Find intervals
intervals = [];
for i = 1:(length(x_vals)-1)
    if f(x_vals(i)) * f(x_vals(i+1)) < 0
        intervals = [intervals; x_vals(i), x_vals(i+1)];
    end
end

% Find roots
roots = [];
for i = 1:size(intervals, 1)
    roots = [roots, hybrid(@f, @dfdx, intervals(i, 1), intervals(i, 2), 1e-12, 1e-16)];
end

roots


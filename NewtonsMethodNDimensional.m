function val = f(x)
    val = [x(1)^2 + x(2)^4 + x(3)^6 - 2;
        cos(x(1) * x(2) * x(3)^2) - x(1) - x(2) - x(3);
        x(2)^2 + x(3)^3 - (x(1) + x(2) - x(3))^2];
end

function j = jac(x)
   j = [ 2*x(1) , 4*x(2)^3, 6*x(3)^5; 
        -x(2)*(x(3)^2)*sin(x(1) * x(2) * x(3)^2) - 1, -x(1)*(x(3)^2)*sin(x(1) * x(2) * x(3)^2) - 1, -2*x(1)*x(2)*x(3)*sin(x(1) * x(2) * x(3)^2) - 1;
        -2*(x(1)+x(2)-x(3)), 2*(x(3)-x(1)), 2*x(1) + 2*x(2) + x(3)*(3*x(3) -2)];
end

x_0 = [ -1.00 ; 0.75 ; 1.50 ];

function x = newtond(f, jac, x0, tol)
% f: Function which implements the nonlinear system of equations.
% Function is of the form f(x) where x is a length-d vector, and
% which returns length-d column vector.
% jac: Function which is of the form jac(x) where x is a length-d vector, and
% which returns the d x d matrix of Jacobian matrix elements for the
% nonlinear system defined by f.
% x0: Initial estimate for iteration (length-d column vector).
% tol: Convergence criterion: routine returns when relative magnitude
% of update from iteration to iteration is <= tol.
% Output => x: Estimate of root (length-d column vector)
    max_iterations = 50;
    iteration = 0;
    converged = false;
    
    dx = jac(x0)\f(x0);
    
    x = x0;

    while ~converged

        if iteration > max_iterations   % Check to see if too many iterations have occured
            warning("Exceeded max iterations, breaking with current values")
            break
        end
        
        % Perform iteration
        dx = jac(x)\f(x);
        x = x - dx;

        % Relative tolerance
        if rms(dx)/rms(x) <= tol
            converged = true;
        end

        iteration = iteration + 1;

    end

end


root = newtond(@f, @jac, x_0, 1e-16)

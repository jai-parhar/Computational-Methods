function yout = rk4step(fcn, t0, dt, y0)
% Performs single iteration of RK4 algorithm to solve ODEs
% Inputs
% fcn: Function handle for right hand sides of ODEs (returns
% length-n column vector).
% t0: Initial value of independent variable.
% dt: Time step.
% y0: Initial values (length-n column vector).
%
% Output
% yout: Final values (length-n column vector)

% Calculate RK4 f values
f0 = fcn(t0, y0);
f1 = fcn(t0 + (dt/2), y0 + (dt/2)*f0);
f2 = fcn(t0 + (dt/2), y0 + (dt/2)*f1);
f3 = fcn(t0 + dt, y0 + dt*f2);

% Perform single RK4 step
yout = y0 + (dt/6)*(f0 + 2*f1 + 2*f2 + f3);

end
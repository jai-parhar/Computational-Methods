# Computational Methods

Various computational methods with applications to mathematical and physical simulations. Written in MATLAB. Primarily written over the course of UBC course PHYS-410, with minor modifications since then. Occasioanlly used in smaller testing simulations since.

## File contents
This will explain the contents of each file at a high level. Further documentation exists within each source file, at the beginning of each function.

#### HybridRootFinder.m
A hybrid method to find roots for 1 dimensional functions. Performs bisection to find location close to a root within tolerance `tol1`, then performs Newton's method within tolerance `tol2` to find the root.

#### NewtonsMethodNDimensional.m
A root finder for N dimensional functions. Performs Newton's method to find roots within tolerance `tol1`.

### RK4 ODE Solver

#### rk4.m
Implementation of a fixed-timestep fourth order Runge-Kutta solver. Works with any order of ODE.

#### rk4ad.m
Implementation of a variable-timestep fourth order Runge-Kutta solver. Automatically picks timesteps to keep error within relative tolerance `reltol`. Works with any order of ODE.

#### rk4step.m
Implementation of a single fourth order Runge-Kutta iteration. Required for use of other files.

#### Applications and Error Analysis
Use of the previous files on specific ODEs, and error analysis code.

##### trk4_sho.m
Application and error analysis of fixed-timestep RK4 method on simple harmonic oscillator ODE (x'' = -x).

##### trk4ad_sho.m
Application and error analysis of variable-timestep RK4 method on simple harmonic oscillator ODE (x'' = -x).
	
##### trk4_vdp.m
Application of fixed-timestep RK4 method on Van der Pol oscillator ODE (x'' + a(x^2 - 1)x' + x = 0).
	
##### trk4ad_vdp.m
Application of variable-timestep RK4 method on Van der Pol oscillator ODE (x'' + a(x^2 - 1)x' + x = 0).

##### trk4step.m
Testing of the error of a single RK4 step on simple harmonic oscillator ODE (x'' = -x).


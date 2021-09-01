# dirk-nonlinear-ode-solver

Solver for nonlinear ODEs of the form
<img src="https://latex.codecogs.com/gif.latex?
	\begin{cases}
		y'(t) = f(y,t) \\
		y(0)  = y_0 
	\end{cases} 
" /> 
using the Diagonally Implicit Runge Kutta methods and Newton solver for the solution of nonlinear equations.

### Usage
- Choose 
    - the desirable function `f`, 
    - its derivative <img src="https://latex.codecogs.com/gif.latex? \partial_y f" /> (`fprime`), 
    - initial condition `yO` and 
    - the exact solution `exact`
- Choose
	- Number of the time step refinements `NT`,
	- final time `T`    
	- tolerance for the Newton solver `tol`

- The output is the convegence plot for each method

### Implemented methods
- DIRK(2,2): 2 stages, 2nd order, [1]
- DIRK(3,3): 3 stages, 3rd order, [2]
- DIRK(5,4): 5 stages, 4th order, [3]

### Tested with functions:
- <img src="https://latex.codecogs.com/gif.latex? f = y+t, \qquad y_0 = 0 " /> 
- <img src="https://latex.codecogs.com/gif.latex? f = yt, \qquad y_0 = 1 " />
- <img src="https://latex.codecogs.com/gif.latex? f = t^n,\qquad n \in \mathbb{N}^+ \qquad y_0 = 1 " />  
- <img src="https://latex.codecogs.com/gif.latex? f = -10000\left( y-\sin(t+\pi/4) \right) + \cos(t + \pi/4), \qquad y_0 = \sin(\pi/4)" /> 


### References

[1] R. Alexander. Diagonally Implicit Runge–Kutta Methods for Stiff O.D.E.’s. _SIAM Journal on Numerical Analysis_, 14(6):1006–1021, 1977.

[2] J. R. Cash. Diagonally implicit Runge-Kutta formulae with error estimates. _IMA Journal of Applied Mathematics (Institute of Mathematics and Its Applications)_,
24(3):293–301, 1979.

[3] G. Wanner and E. Hairer. _Solving Ordinary Differential Equations II. Stiff and Differential-Algebraic Problems._ Springer-Verlag, Berlin Heidelberg, 1996.
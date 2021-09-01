# dirk-nonlinear-ode-solver

Solver for nonlinear ODEs of the form

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{cases}&space;y'(t)&space;=&space;f(y,t)&space;\\&space;y(0)&space;=&space;y_0&space;\end{cases}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{cases}&space;y'(t)&space;=&space;f(y,t)&space;\\&space;y(0)&space;=&space;y_0&space;\end{cases}" title="\begin{cases} y'(t) = f(y,t) \\ y(0) = y_0 \end{cases}" /></a>

using the Diagonally Implicit Runge Kutta methods and Newton solver for the solution of nonlinear equations.

### Usage
- Choose 
    - the desirable function `f`, 
    - its derivative <a href="https://www.codecogs.com/eqnedit.php?latex=\partial_y&space;f" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\partial_y&space;f" title="\partial_y f" /></a> (`fprime`), 
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
- <a href="https://www.codecogs.com/eqnedit.php?latex=f&space;=&space;y&plus;t,&space;\qquad&space;y_0&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f&space;=&space;y&plus;t,&space;\qquad&space;y_0&space;=&space;0" title="f = y+t, \qquad y_0 = 0" /></a> 
- <a href="https://www.codecogs.com/eqnedit.php?latex=f&space;=&space;yt,&space;\qquad&space;y_0&space;=&space;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f&space;=&space;yt,&space;\qquad&space;y_0&space;=&space;1" title="f = yt, \qquad y_0 = 1" /></a>
- <a href="https://www.codecogs.com/eqnedit.php?latex=f&space;=&space;t^n,\qquad&space;n&space;\in&space;\mathbb{N}^&plus;&space;\qquad&space;y_0&space;=&space;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f&space;=&space;t^n,\qquad&space;n&space;\in&space;\mathbb{N}^&plus;&space;\qquad&space;y_0&space;=&space;1" title="f = t^n,\qquad n \in \mathbb{N}^+ \qquad y_0 = 1" /></a> 
- <a href="https://www.codecogs.com/eqnedit.php?latex=f&space;=&space;-10000\left(&space;y-\sin(t&plus;\pi/4)&space;\right)&space;&plus;&space;\cos(t&space;&plus;&space;\pi/4),&space;\qquad&space;y_0&space;=&space;\sin(\pi/4)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f&space;=&space;-10000\left(&space;y-\sin(t&plus;\pi/4)&space;\right)&space;&plus;&space;\cos(t&space;&plus;&space;\pi/4),&space;\qquad&space;y_0&space;=&space;\sin(\pi/4)" title="f = -10000\left( y-\sin(t+\pi/4) \right) + \cos(t + \pi/4), \qquad y_0 = \sin(\pi/4)" /></a>


### References

[1] R. Alexander. Diagonally Implicit Runge–Kutta Methods for Stiff O.D.E.’s. _SIAM Journal on Numerical Analysis_, 14(6):1006–1021, 1977.

[2] J. R. Cash. Diagonally implicit Runge-Kutta formulae with error estimates. _IMA Journal of Applied Mathematics (Institute of Mathematics and Its Applications)_,
24(3):293–301, 1979.

[3] G. Wanner and E. Hairer. _Solving Ordinary Differential Equations II. Stiff and Differential-Algebraic Problems._ Springer-Verlag, Berlin Heidelberg, 1996.
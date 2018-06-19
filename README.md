# Averaged Control

In this work, we address the optimal control of parameter-dependent systems. We introduce the notion of averaged control in which the quantity of interest is the average of the states with respect to the parameter family <img src="https://latex.codecogs.com/gif.latex?%5Cmathcal%7BK%7D%3D%20%5Cleft%5C%7B%20%5Cnu_i%20%5Cin%20%5Cmathbb%7BR%7D%2C%20%5Censpace%201%5Cleq%20i%20%5Cleq%20K%20%5Cright%5C%7D">.

<p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cmin%20_%7Bu%20%5Cin%20%5Cmathcal%7BU%7D_%7Bad%7D%7D%20%5Cmathcal%7BJ%7D%5Cleft%28%20u%5Cright%29%20%3D%20%5Cmin%20_%7Bu%20%5Cin%20%5Cmathcal%7BU%7D_%7Bad%7D%7D%20%5Cfrac%7B1%7D%7B2%7D%20%5Cleft%5B%20%5Cfrac%7B1%7D%7BK%7D%20%5Csum_%7B%5Cnu%20%5Cin%20%5Cmathcal%7BK%7D%7D%20x%20%5Cleft%28%20T%2C%20%5Cnu%20%5Cright%29%20-%20%5Cbar%7Bx%7D%20%5Cright%5D%5E2%20&plus;%20%5Cfrac%7B%5Cbeta%7D%7B2%7D%20%5Cint_0%5ET%20u%5E2%20%5Cmathrm%7Bd%7Dt%2C%20%5Cquad%20%5Cbeta%20%5Cin%20%5Cmathbb%7BR%7D%5E&plus;%2C">
</p>

with <img src="https://latex.codecogs.com/gif.latex?%5Cmathcal%7BU%7D_%7Bad%7D"> the space of admissible controls and <img src="https://latex.codecogs.com/gif.latex?%5Cbar%7Bx%7D"> the average state target. The optimization problem is subject to the finite dimensional linear control system

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20x%5E%5Cprime%20%5Cleft%28%20t%20%5Cright%29%20%3D%20A%20%5Cleft%28%20%5Cnu%20%5Cright%29%20x%20%5Cleft%28%20t%20%5Cright%29%20&plus;%20B%20%5Cleft%28%20%5Cnu%20%5Cright%29%20u%20%5Cleft%28%20t%20%5Cright%29%2C%20%5Cquad%200%20%3C%20t%20%3CT%2C%20%5C%5C%20x%7B%5Cleft%28%200%20%5Cright%29%7D%20%3D%20x%5E0.%20%5Cend%7Bcases%7D">
</p>

## Steepest Descent Method

We use the classical gradient descent method based on the adjoint methodology, and obtain the corresponding adjoint system,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20p%5E%5Cprime%20%5Cleft%28%20t%20%5Cright%29%20%3D%20-A%20%5Cleft%28%20%5Cnu%20%5Cright%29%20p%20%5Cleft%28%20t%20%5Cright%29%2C%20%5Cquad%200%20%3C%20t%20%3CT%2C%20%5C%5C%20p%7B%5Cleft%28%20T%20%5Cright%29%7D%20%3D%20-%20%5Cleft%5B%20%5Cdisplaystyle%20%5Cfrac%7B1%7D%7BK%7D%20%5Csum_%7B%5Cnu%20%5Cin%20%5Cmathcal%7BK%7D%7D%20x%20%5Cleft%28%20T%2C%20%5Cnu%20%5Cright%29%20-%20%5Cbar%7Bx%7D%20%5Cright%5D.%20%5Cend%7Bcases%7D">
</p>

The functional is minimized by taking the steepest descent direction given by

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?J%5E%5Cprime%20%5Cleft%5B%20u%5E%7B%5Cleft%28%20k%20%5Cright%29%7D%20%5Cright%5D%20%3D%20%5Cbeta%20u%5E%7B%5Cleft%28%20k%20%5Cright%29%7D%5Cleft%28%20t%20%5Cright%29%20-%20%5Cfrac%7B1%7D%7BK%7D%20%5Csum_%7B%5Cnu%20%5Cin%20%5Cmathcal%7BK%7D%7D%20B%5Et%20%5Cleft%28%20%5Cnu%20%5Cright%29%20p%5E%7B%5Cleft%28%20k%20%5Cright%29%7D%20%5Cleft%28%20t%2C%20%5Cnu%20%5Cright%29,">
</p>

and the new control reads as

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?u%5E%7B%5Cleft%28%20k&plus;1%20%5Cright%29%7D%20%3D%20u%5E%7B%5Cleft%28%20k%20%5Cright%29%7D%20-%20%5Cepsilon%20J%5E%5Cprime%20%5Cleft%5B%20u%5E%7B%5Cleft%28%20k%20%5Cright%29%7D%5Cright%5D">
</p>

for some <img src="https://latex.codecogs.com/gif.latex?%5Cepsilon"> small enough.

## Conjugate Gradient Method

We have also used the conjugate gradient method in order to reach faster the optimal control. In order to be able to apply this method the state vector has been split as

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?x%5Cleft%28%20t%20%5Cright%29%20%3D%20z_u%20%5Cleft%28%20t%20%5Cright%29%20&plus;%20y%5Cleft%28%20t%20%5Cright%29%2C">
</p>

where <img src="https://latex.codecogs.com/gif.latex?z_u"> is the solution to the controlled system with zero initial condition,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20z%5E%5Cprime_u%20%5Cleft%28%20t%20%5Cright%29%20%3D%20A%20%5Cleft%28%20%5Cnu%20%5Cright%29%20z_u%20%5Cleft%28%20t%20%5Cright%29%20&plus;%20B%20%5Cleft%28%20%5Cnu%20%5Cright%29%20u%20%5Cleft%28%20t%20%5Cright%29%2C%20%5Cquad%200%20%3C%20t%20%3CT%2C%20%5C%5C%20z_u%7B%5Cleft%28%200%20%5Cright%29%7D%20%3D%200%2C%20%5Cend%7Bcases%7D">
</p>

and <img src="https://latex.codecogs.com/gif.latex?y"> solves the free dynamics problem,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20y%5E%5Cprime%20%5Cleft%28%20t%20%5Cright%29%20%3D%20A%20%5Cleft%28%20%5Cnu%20%5Cright%29%20y%20%5Cleft%28%20t%20%5Cright%29%2C%20%5Cquad%200%20%3C%20t%20%3CT%2C%20%5C%5C%20y%7B%5Cleft%28%200%20%5Cright%29%7D%20%3D%20x%5E0.%20%5Cend%7Bcases%7D">
</p>

The functional can be expressed as

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?J%5Cleft%28%20u%20%5Cright%29%20%3D%20%5Cfrac%7B1%7D%7B2%7D%5Cleft%28%20%5Cbar%7Bz%7D_u%5Cleft%28T%5Cright%29%20&plus;%20%5Cbar%7By%7D%20%5Cleft%28%20T%20%5Cright%29%20-%20%5Cbar%7Bx%7D%2C%20%5Cbar%7Bz%7D_u%5Cleft%28T%5Cright%29%20&plus;%20%5Cbar%7By%7D%20%5Cleft%28%20T%20%5Cright%29%20-%20%5Cbar%7Bx%7D%20%5Cright%29_%7B%5Cdisplaystyle%20%5Cmathbb%7BR%7D%5En%7D%20&plus;%20%5Cfrac%7B%5Cbeta%7D%7B2%7D%20%5Cleft%28u%2C%20u%5Cright%29_%7BL%5E2%5Cleft%28%5Cleft%5B0%2CT%5Cright%5D%5Cright%29%7D.">
</p>

We introduce the linear operator

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%5CLambda%3A%20L%5E2%20%5Cleft%28%20%5B0%2C%20T%5D%20%5Cright%29%20%26%20%5Crightarrow%20%5Cmathbb%7BR%7D%5En%20%5C%5C%20u%20%26%20%5Crightarrow%20%5Cbar%7Bz%7D_u%20%5Cleft%28%20T%20%5Cright%29%20%3D%20%5Cfrac%7B1%7D%7BK%7D%20%5Csum_%7B%5Cnu%20%5Cin%20%5Cmathcal%7BK%7D%7D%20z_u%20%5Cleft%28%20T%2C%20%5Cnu%20%5Cright%29%20%5Cend%7Balign*%7D">
</p>

and its dual counterpart,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%5CLambda%5E*%3A%20%5Cmathbb%7BR%7D%5En%20%26%20%5Crightarrow%20L%5E2%20%5Cleft%28%20%5B0%2C%20T%5D%20%5Cright%29%20%5C%5C%20p_T%20%26%20%5Crightarrow%20-%20%5Cfrac%7B1%7D%7BK%7D%20%5Csum_%7B%5Cnu%20%5Cin%20%5Cmathcal%7BK%7D%7D%20B%5Et%20%5Cleft%28%20%5Cnu%20%5Cright%29%20p%20%5Cleft%28%20t%2C%20%5Cnu%20%5Cright%29%20%5Cend%7Balign*%7D">
</p>

where <img src="https://latex.codecogs.com/gif.latex?p%5Cleft%28t%2C%20%5Cnu%5Cright%29"> is solution to

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20p%5E%5Cprime%20%5Cleft%28%20t%20%5Cright%29%20%3D%20-A%20%5Cleft%28%20%5Cnu%20%5Cright%29%20p%20%5Cleft%28%20t%20%5Cright%29%2C%20%5Cquad%200%20%3C%20t%20%3CT%2C%20%5C%5C%20p%7B%5Cleft%28%20T%20%5Cright%29%7D%20%3D%20-%20p_T.%20%5Cend%7Bcases%7D">
</p>

By doing this we can write the directional derivative of the functional as

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cmathcal%7BD%7D_%7B%5Cdelta%20u%20%7DJ%20%5Cleft%28%20u%20%5Cright%29%20%3D%20%5Cleft%28%20%5Cunderbrace%7B%5Cleft%28%20%5CLambda%5E*%20%5CLambda%20&plus;%20%5Cbeta%20I%20%5Cright%29%7D_%7B%3DA_%7Bcg%7D%7D%20u%20-%20%5Cunderbrace%7B%5CLambda%5E*%5Cleft%28%20%5Cbar%7Bx%7D%20-%20%5Cbar%7By%7D%20%5Cleft%28%20T%20%5Cright%29%20%5Cright%29%7D_%7B%3Db_%7Bcg%7D%7D%2C%20%5Cdelta%20u%20%5Cright%29_%7BL%5E2%5Cleft%28%5Cleft%5B0%2CT%5Cright%5D%5Cright%29%7D.">
</p>

After having defined <img src="https://latex.codecogs.com/gif.latex?A_%7Bcg%7D"> and <img src="https://latex.codecogs.com/gif.latex?b_%7Bcg%7D"> we can apply the conjugate gradient method to solve the control problem.


## Running the example

Both steepest descent and conjugate gradient methods have been implemented in Matlab. The can be run by typing in the Command Window

* for the steepest descent method:
```Matlab
AveragedControlSD
```

* for the conjugate gradient method:
```Matlab
AveragedControlCG
```

<p align="center">
  <img src="u.png">
</p>

<p align="center">
  <img src="xav.png">
</p>

<p align="center">
  <img src="xi.png">
</p>

## References

* E. Zuazua (2014). Averaged Control. _Automatica_, 50 (12), p. 3077-3087.

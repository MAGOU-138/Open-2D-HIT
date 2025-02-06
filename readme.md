# Direct numerical simulation of 2-dimensional homogeneous isotropic turbulence using Fourier spectral methods
![][shot-198.png][]

## keywords
- spectral methods, Fourier, DNS, Kraichnan turbulence, homogeneous isotropic turbulence (HIT)

# Features
- OpenMPI parallelization, one direction decomposition
- FFTW used for Fourier modes
- High-order time integration scheme
- Lighweight, low encapsulation level, easy to read and modify

# Governing equations

$$
\begin{aligned}
\partial_t\omega+\frac{\partial \psi}{\partial y}\frac{\partial \omega}{\partial x}-\frac{\partial \psi}{\partial x}\frac{\partial \omega}{\partial y}=\nu\nabla^2\omega-\alpha\omega+f,\\
\omega=-\nabla^2\psi.
\end{aligned}
$$

This equation solves the scalar vorticity $\omega$ and streamfunction $\psi$. The $\alpha$ is a linear frictional damping, and $f$ is the curl of forcing term.

<!-- The dimension of vorticity, streamfunction, and other parameters are listed in table. -->

<!-- |$\omega$|$\psi$|$\nu$|$\alpha$|$f$|
|:---:|:---:|:---:|:---:|:---:|
|$[T^{-1}]$|$[L^2T^{-1}]$|$[L^2T^{-1}]$|$[T^{-1}]$|$[T^{-2}]$| -->

<!-- The flow is confined in a cyclic box of side $L \times  L$. We expand the vorticty and streamfunction in Fourier series so that the equation becomes -->

$$
\omega(x_i,y_i)=\sum_{m=-N_x/2}^{N_x/2-1}\sum_{n=-N_y/2}^{N_y/2-1}\hat{\omega}_{m,n}\exp\left[\mathrm{i}\left(\frac{2\pi m}{L_x}x_i+\frac{2\pi n}{L_y}y_j\right)\right].
$$

<!-- We note

$$
k_x=\frac{2\pi m}{L_x},\,k_y=\frac{2\pi n}{L_y}.
$$ -->

In Fourier wavenumber space, the equation becomes

$$
\left[\frac{\partial}{\partial t}+\nu k^2+\alpha\right]\hat{\omega}_{m,n}=\hat{J},
$$

where $J=\frac{\partial \psi}{\partial x}\frac{\partial \omega}{\partial y}-\frac{\partial \psi}{\partial y}\frac{\partial \omega}{\partial x}$ is the non-linear term.

The time integration scheme is a multi-step scheme [1]

$$
\frac{\partial\psi^{(n+1)}}{\partial t}\Delta t=\frac{11}{6}\psi^{(n+1)}-\left(3\psi^{(n)}-\frac{3}{2}\psi^{(n-1)}+\frac{1}{3}\psi^{(n-2)}\right),
$$

and the nonlinear term is approximated by Taylor expansion

$$
\hat{J}^{(n+1)}=3\hat{J}^{(n)}-3\hat{J}^{(n-1)}+\hat{J}^{(n-2)},
$$

and the viscous terms are treated implicitly.


The discrete algebra equation is then

$$
\left[\gamma_0+\Delta t\left(\nu k^2+\alpha\right)\right]\hat{\omega}^{(n+1)}=\sum_{q=0}^{2}\alpha_q\omega^{(n-q)}+\Delta t\sum_{q=0}^{2}\beta_q\hat{J}^{(n-q)},
$$

where
|$\gamma_0$|$\alpha_0$|$\alpha_1$|$\alpha_2$|$\beta_0$|$\beta_1$|$\beta_2$|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|$11/6$|$3$|$-3/2$|$1/3$|$3$|$-3$|$1$|

# Installation
## Prerequisites
- MPI
- HDF5
- WFFT


## Complilation
Modify the `makefile` with your own libraries' pathes, and then `make`, just as easy as that.


# Reference
- [1] G. E. Karniadakis, M. Israeli, S. A. Orszag, High-order splitting methods for the incompressible Navier-Stokes equations, Journal of Computational Physics 97 (2) (1991) 414–443. doi:https://doi.org/12410.1016/0021-9991(91)90007-8.


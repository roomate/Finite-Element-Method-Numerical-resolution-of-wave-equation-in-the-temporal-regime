## Numerical resolution of wave equation in the temporal regime

This mini-project is a standard application of Finite Element Method for elliptic and hyperbolic PDEs with Dirichlet boundary condition. The first part concerns the resolution of Laplace equation with a non-homogeneous elasticity in 
a rectangular domain while in a second time, one deals with the propagation of a wave in this same inhomogeneous medium.

### Stationary resolution
The resolution of the Laplace equation via the Finite Element Method is standard here. After an adequate discretization of the Sobolev space, it amounts to solving a linear equation. The associated matrix is known to be sparse; MATLAB takes advantage of such a structure to accelerate the resolution of the linear equation. 

<img src="img/Solexact.jpg" alt="drawing" width="400"/>

### Temporal resolution
The temporal discretization is made with a leap-frog scheme; three different techniques are employed for the space discretization: 
- Standard Quadrature,
- Mass condensation,
- and Cholesky decoomposition.

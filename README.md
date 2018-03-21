
**Author:** James Rekow
**Overview:** This code contains an example of a boundary value problem (differential equation with boundary conditions) and an interface in which the maximum principle is violated but a comparison principle still holds. This is surprising because the traditional proof of the comparison principle for standard boundary value problems uses the maximum principle.

Pure diffusion equation (in u) on (0,1) with Dirichlet BC and an interface
at 0.5. The initial condition is u(x,0) = 1, if x < .5, u(x,0) = 0, x > 0.5.
The interface conditions are (u_{-})_{x}(0.5-) = 0.25(u_{+})_{x}(0.5+)
and u_{-}(0.5-) = 0.25u_{+}(0.5+).

Consider v = u if x < .5, v = 0.25u if x >0.5.
Then equivalent problem in v is a pure diffusion with
Dirichlet BC, no interface, and the same initial data.

In our notation this pair of problems corresponds to 
eta_{-} = kappa_{-} = 1, eta_{+} = kappa_{+} = 0.25.


Note that the interface problem in u must obey a comparison principle,
since the problem in v does, and multiplication by piecewise constants
does not change the ordering of solutions. However, the numerical example
clearly shows that u does not obey a maximum principle. The yellow star in the
time evolving solution to the u problem represents the maximum value of
the solution. It does not attain this maximum on the parabolic boundary.
This is also evidenced by the surface in figure2. The ridge is the value
of u(0.5+,t).

Function returns the point (x,t) at which the solution u to the interface
problem attains its maximum. The red square on the surface marks the
local maximum of u.

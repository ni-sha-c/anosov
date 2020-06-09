using LinearAlgebra
"""
Smoothly perturbed Sawtooth map
Input: 
    x1: array size mx1, m initial conditions
    s: float, norm of perturbation
    n: int, number of timesteps
Output:
    x: array size (n+1)xm
"""
function step(x1, s=0., n=1)
    m, d = size(x1)
    x = zeros(m,n+1) 
    x[1:m] = x1
    f(x) = (2*x + s*sin(2*pi*x)) % 1
    for i=1:n
        x[:,i+1] = f.(x[:,i]) 
    end
    return x'
end
"""
Construct two matrices Phi and Phi1 := Phi circ varphi
of size mxn with elements 
Phi[j,i] = f circ varphi^{i-1}(x[j,:])
Phi1[j,i] (by definition) = f circ varphi^{i}(x[j,:])  

Inputs: 
    x: mxd m-length array of initial conditions each of 
        dimension d
    f: a scalar function which takes x[i] and returns a 
        complex number
    n: number of rows of Phi or time delay
    s: system parameters to be passed to varphi
"""
function construct_hankel_basis(x, s, n, f)
    m, d = size(x)
    x_trj = step(x, s, n)'
    f_trj = f.(x_trj)
    Phi = view(f_trj, 1:m, 1:n)
    Phi1 = view(f_trj, 1:m, 2:(n+1))
    return Phi, Phi1
end
function fun(x)
    return exp(2*pi*im*x)
end
"""
Solves the following minimization problem for H:
    H = arg min_{C in C^{nxn}} ||Phi C - Phi1||^2.

where Phi is an mxn matrix containing evaluations 
of f, f circ varphi, ..., f circ varphi^{n-1} 
at x[1],x[2],..., x[m], and Phi1 := Phi circ varphi.

The solution here is obtained directly as:
         H = (Phi' Phi)^{-1} Phi' Phi1, 
after Phi and Phi' are computed by `construct_hankel_basis`.

Inputs:
    f: a complex-valued function of x[i]
    x: mxd m-length array of d-dimensional states
    n: number of basis functions, Int64
    s: system parameters

"""
function dmd(x, s, n, f)
    if length(f) == 1
        Phi, Phi1 = construct_hankel_basis(x, s, 
    										   n, f[1])
    else
        Phi, Phi1 = construct_basis(x, s, n, f)
    end
    return (Phi' * Phi)\(Phi' * Phi1)
end
"""
Inputs: 
    spe: must be one of the spaces supported by ApproxFun
    n: number of basis functions
"""
function get_basis(spe, n)
    n = 10

end
"""
Plot eigenvalues and eigenvectors of 
H-DMD approximation H of the Koopman operator.

Input: 
    x: mxd matrix of initial conditions each of dimension d
    step: function that takes in x, s(parameters), and n
    outputs varphi^n(x). Output size must be (n+1)xmxd
    

"""
#=
function
=#
  

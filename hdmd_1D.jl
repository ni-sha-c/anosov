include("plucked_map.jl")
using LinearAlgebra
using PyPlot
using ToeplitzMatrices
"""
    Plots eigenvalues of an nxn matrix K representing 
    the Koopman operator

"""
function hdmd()
    n = 15
	n_samples = 2*n
	#x = LinRange(0.,2,n_samples+2)[2:end-1]
	s = 0.01
	m = 0
	x = zeros(n_samples)
	x[1] = 2*rand()
	for i = 1:100
		x[1] = osc_tent(x[1],s,m)
	end
	for i = 2:n_samples
		x[i] = osc_tent(x[i-1],s,m)
	end

	k = 1
	f = exp.(im*pi*k*x)
	F0 = Hankel(f[1:n], f[n:2*n-1])
	F1 = Hankel(f[2:n+1], f[n+1:2*n])
	
	K = (F1*F0')*inv(F0*F0') 
	L = eigvals(K)	
	L_real = real(L)
	L_imag = imag(L)
	return L_real, L_imag, K	
end

L_real, L_imag, K = hdmd()
fig, ax = subplots(1,1)
ax.plot(L_real,L_imag,".")
x = LinRange(0.,2*pi,1000)
ax.plot(cos.(x), sin.(x), "k.", ms=0.1)
ax.plot(0.5*cos.(x), 0.5*sin.(x), "k.", ms=0.1)
ax.plot(0.25*cos.(x), 0.25*sin.(x), "k.", ms=0.1)
ax.plot(0.75*cos.(x), 0.75*sin.(x), "k.", ms=0.1)

#u = rand(1,2)
#s = [0.7,0.3]
#n = 1000


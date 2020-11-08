include("plucked_map.jl")
using LinearAlgebra
using PyPlot
"""
    Plots eigenvalues of an nxn matrix K representing 
    the Koopman operator

"""
function edmd()
    n_basis = 64
	n_samples = 10000
	#x = LinRange(0.,2,n_samples+2)[2:end-1]
	s = 0.5
	m = 10
	n_steps = 500
	x = zeros(n_samples+1)
	x[1] = 2*rand()
	for i = 2:n_samples+1
		x[i] = osc_tent(x[i-1],s,m)
	end
	x1 = x[2:n_samples+1]
	x = x[1:n_samples]
	
	F0 = Array{Complex{Float64}, 2}(undef, n_samples, 
					n_basis) 
	F1 = Array{Complex{Float64}, 2}(undef, n_samples, 
					n_basis)
	K = Array{Complex{Float64}, 2}(undef, n_basis,
					n_basis)
	F0F0 = Array{Complex{Float64}, 2}(undef, n_basis, 
					n_basis) 


	for I in CartesianIndices(F0)
		n, k = Tuple(I)
		F0[I] = exp(pi*im*(k-1)*x[n])
		F1[I] = exp(pi*im*(k-1)*x1[n])
	end
	#=
	for j = 1:n_basis
		F1j = view(F1,:,j)
		F0j = view(F0,:,j)
		for i = 1:n_basis
			F0i = view(F0,:,i)
			K[i,j] = dot(F0i,F1j)/n_samples
			F0F0[i,j] = dot(F0i,F0j)/n_samples 
		end
	end
	K = F0F0'\K'
	=#
	K = (F1'*F0)*inv(F0'*F0) 
	L = eigvals(K)	
	L_real = real(L)
	L_imag = imag(L)
	return L_real, L_imag	
end

L_real, L_imag = edmd()
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


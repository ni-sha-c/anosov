include("plucked_map.jl")
include("wave_map.jl")
using LinearAlgebra
using PyPlot
"""
    Plots eigenvalues of an nxn matrix K representing 
    the Koopman operator

"""
function edmd_wave(s,n_basis)
	n_samples = 10000
	#x = LinRange(0.,2,n_samples+2)[2:end-1]
	n_steps = 500
	x = zeros(n_samples+1)
	x[1] = 2*rand()
	for i = 1:n_steps
		x[1] = wave_tent(x[1],s)
	end
	for i = 2:n_samples+1
		x[i] = wave_tent(x[i-1],s)
	end
	x1 = x[2:n_samples+1]
	x = x[1:n_samples]
	
	F0 = Array{Complex{Float64}, 2}(undef, n_samples, 
					n_basis) 
	F1 = Array{Complex{Float64}, 2}(undef, n_samples, 
					n_basis)
	for I in CartesianIndices(F0)
		n, k = Tuple(I)
		F0[I] = exp(pi*im*(k-1)*x[n])
		F1[I] = exp(pi*im*(k-1)*x1[n])
	end
	return solveK(F0,F1)
end
function edmd_pluck(s,m,n_basis)
	n_samples = 10000
	n_steps = 500
	x = zeros(n_samples+1)
	x[1] = 2*rand()
	for i = 1:n_steps
		x[1] = osc_tent(x[1],s,m)
	end
	for i = 2:n_samples+1
		x[i] = osc_tent(x[i-1],s,m)
	end
	x1 = x[2:n_samples+1]
	x = x[1:n_samples]
	
	F0 = Array{Complex{Float64}, 2}(undef, n_samples, 
					n_basis) 
	F1 = Array{Complex{Float64}, 2}(undef, n_samples, 
					n_basis)

	for I in CartesianIndices(F0)
		n, k = Tuple(I)
		F0[I] = exp(pi*im*(k-1)*x[n])
		F1[I] = exp(pi*im*(k-1)*x1[n])
	end
	return solveK(F0, F1)
end


function solveK(F0,F1)
	K =	 (F1'*F0)*inv(F0'*F0) 
	L = eigvals(K)	
	L_real = real(L)
	L_imag = imag(L)
	return L_real, L_imag
end



s = 0.001 
m = 0
n = 64
L_real, L_imag = edmd_pluck(s, m, n)
fig, ax = subplots(1,1)
ax.plot(L_real,L_imag,"o",label="s = $s, n = $m")


s = 0.1
m = 6
L_real, L_imag = edmd_pluck(s, m, n)
ax.plot(L_real,L_imag,"<",label="s = $s, n = $m")


x = LinRange(0.,2*pi,1000)
ax.plot(cos.(x), sin.(x), "k.", ms=0.1)
ax.plot(0.5*cos.(x), 0.5*sin.(x), "k.", ms=0.1)

ax.legend(fontsize=20)
ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)
ax.grid(true)
ax.axis("scaled")
#u = rand(1,2)
#s = [0.7,0.3]
#n = 1000


include("step.jl")
using LinearAlgebra
using PyPlot
"""
    Plots eigenvalues of an nxn matrix K representing 
    the Koopman operator

"""
#function edmd()
    n_basis = 5
	n_basis_tot = n_basis*n_basis
	n_grid = 10
	x = LinRange(0.,1,n_grid+1)[1:end-1]
	n_samples = n_grid*n_grid
	x0_g = repeat(x, n_grid, 1)
	y0_g = reshape(repeat(x',n_grid,1),n_samples,:)
	u0_g = reshape([x0_g; y0_g], n_samples, 2) 
	s = [0.7,0.3]
	u_g = step(u0_g, s, 1, n_samples)
	u_g_0 = view(u_g, :, :, 1)'
	u_g_1 = view(u_g, :, :, 2)'
	F0 = Array{Complex{Float64}, 3}(undef, n_samples, 
					n_basis, n_basis) 
	F1 = Array{Complex{Float64}, 3}(undef, n_samples, 
					n_basis, n_basis)
	K = Array{Complex{Float64}, 2}(undef, n_basis_tot,
					n_basis_tot)

	for I in CartesianIndices(F0)
		n, k1, k2 = Tuple(I)
		F0[I] = exp(2*pi*im*dot([k1,k2],u_g_0[:,n]))
		F1[I] = exp(2*pi*im*dot([k1,k2],u_g_1[:,n]))
	end
	Ind = CartesianIndices((1:n_basis,1:n_basis))
	for J in Ind
		F1j = view(F1,:,J)
		j = LinearIndices(Ind)[J] 
		for I in Ind
			F0i = view(F0,:,I)
			i = LinearIndices(Ind)[I]
			K[i, j] = dot(F0i,F1j)/n_samples
		end
	end
	L = eigvals(K)	
	L_real = real(L)
	L_imag = imag(L)
	fig, ax = subplots(1,1)
	ax.plot(L_real,L_imag,".")
	#ax.plot(x, sqrt.(1.0 .- x.*x), "k.")


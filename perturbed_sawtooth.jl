"""
Smoothly perturbed Sawtooth map
Input: 
	x1: array size m, m initial conditions
	s: float, norm of perturbation
	n: int, number of timesteps
Output:
	x: array size nxm
"""
function step(x1, s=0., n=1)
	m, = size(x1)
	x = zeros(m,n) 
	x[1:m] = x1
	f(x) = (2*x + s*sin(2*pi*x)) % 1
	for i=1:n
		@time x[:,i] = f.(x[:,i-1]) 
	end
	return x'
end
"""
Hankel DMD:
	
"""
function hdmd(x, f)
	n, = size(x)
	inds = axes(x,1)
	fx = zeros(n,n)
	fstepx = zeros(n,n)
	for i = inds 
		fx[:,i] = f.(x)
	end
	
		

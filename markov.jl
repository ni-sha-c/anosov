function ulam(x,nbins)
	P = zeros(nbins, nbins)
	dx = 2.0/nbins
	n = length(x)
	for k = 1:n-1
		binc = ceil(Int64,x[k]/dx)
		binn = ceil(Int64,x[k+1]/dx)
		P[binc,binn] += 1/n
	end
	return P
end



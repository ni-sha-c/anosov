function markov(x,nbins)
	P = zeros(nbins, nbins)
	dx = 2.0/nbins
	n = length(x)
	for k = 1:n-1
		binc = floor(Int64,x[k]/dx)+1
		binn = floor(Int64,x[k+1]/dx)+1
		P[binc,binn] += 1/n
	end
	for i = 1:nbins
		sp = sum(P[i,:])
		if sp != 0
			P[i,:] ./= sp
		end
	end
	return P
end



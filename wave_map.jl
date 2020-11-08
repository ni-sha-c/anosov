function wave_tent(x,s)
    if x < 1
		pert = sqrt( (1+s)^2 - 4*s*x)
		return 4*x/(1 + s + pert)
	end
	pert = sqrt( (1-s)^2 - 4*s*(2-x))
	return  4*(2-x)/(1 - s + pert)
end


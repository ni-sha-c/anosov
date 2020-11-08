function tent_basic(x,s)
    if x < 1
        return min(2*x/(1-s), 
                2 - 2*(1-x)/(1+s))
	end
    return min(2*(2-x)/(1-s), 
            2 - 2*(x-1)/(1+s))
end
function oscillation(x,s)
    if x < 0.5
        return tent_basic(2*x,s)/2
	end
    return 2-tent_basic(2-2*x,s)/2
end
function frequency(x,s,n)
    return oscillation(2^n*x - floor(2^n*x),s)/2^n + 
            2*floor(2^n*x)/2^n
end
function osc_tent(x, s, n)
    return min(frequency(x,s,n), frequency(2-x,s,n))
end


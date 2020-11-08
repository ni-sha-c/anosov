using PyPlot
include("plucked_map.jl")
s = [0.0, 0.1, 0.2, 0.5]
n = [0, 3, 6]
x = LinRange(0.,2.0,100)
for nk in n
	fig = figure(figsize=(8,6))
	ax = fig.add_subplot(111)
    for sk in s
    	ax.plot(x, osc_tent.(x, sk, nk), label="n = $(nk), s = $(sk)")
	end
	ax.legend(fontsize=18)
	ax.xaxis.set_tick_params(labelsize=18)
	ax.yaxis.set_tick_params(labelsize=18)
end




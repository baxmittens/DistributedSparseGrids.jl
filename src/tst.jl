include("AdaptiveSparseGrids.jl")

N = 1
CT = Float64
RT = Float64
CPType = CollocationPoint{N,CT}

HCPType = HierarchicalCollocationPoint{N,CPType,RT}

Maxp = 1
nvmom=Maxp
Nv = numwaveletcoeffs(N, nvmom)
force_zero_mean = true
WCPType = WaveletCollocationPoint{N,CPType,RT,CT,Nv}

pointprobs = SVector{N,Int}([1 for i = 1:N])

@time begin
wasg = init(AHSG{N,WCPType},pointprobs,Maxp)
for i = 1:3; generate_next_level!(wasg); end
end
#f = x-> 1/(abs(sqrt(2)-(x[1]-1)^2)+.5)
#f = x-> x[1]< -0.5 ?  x[1]^2.0 : x[1] < 0.1 ? sin(2*x[1]-5) : x[1] < .5 ? 0.01/x[1] : log(x[1])
#f = x-> x[1]^2*x[2]^2
f = x->abs(norm(x)-0.75) * ifelse(norm(x) > 0.75, 0.5, 5.0) * x[2] + pi/100.
#f =  x-> sqrt(sum(x.^2)) > .5 ? 1.0 : -1.0

wavelet_coeffs!(wasg, nvmom)
init_weights!(wasg, f)
wavelet_weights!(wasg,f,nvmom)


wavlet(hcpt,x) = begin
	#nvmom = numlevelsdown_possible(hcpt, numlevelsdown)
	#numlevelsdown = min(numlevelsdown,nvmom)
	res = basis_fun(hcpt, [x], Maxp)
	untillvl = level(hcpt)-nvmom
	aci = AncestorIterator(hcpt, untillvl)
	icpt = 0
	for anc in aci
		icpt += 1
		res -= hcpt.coeffs[icpt] * basis_fun(anc, [x], Maxp)
	end
	return res
end

import PlotlyJS
plots = PlotlyJS.GenericTrace{Dict{Symbol,Any}}[]
xx = linspace(-1,1,100)
for hcpt in wasg
	push!(PlotlyJS.scatter(x=xx,y=map(x->wavlet(hcpt,x),xx)))
end


#for i = 1:3
#	nchilds = generate_next_level!(wasg,1e-3)
#	@info "$(length(nchilds)) new CollocationPoints created"
#	wavelet_coeffs!(wasg, nvmom)
#	init_weights!(wasg, f)
#	wavelet_weights!(wasg,f,nvmom)
#end

#for i = 1:5
#hcpts = wavelet_refine!(wasg, 1e-2)
#println("$(length(hcpts)) new CollocationPoints")
#wavelet_coeffs!(wasg, Maxp-1)
#wavelet_weights!(wasg,f,Maxp-1)
#init_weights!(wasg, f)
#end



display(PlotlyJS.plot([PlotlyJS.scatter3d(wasg,true,Maxp,mode="markers"),PlotlyJS.surface(wasg,200)]))
#display(PlotlyJS.plot([PlotlyJS.scatter3d(wasg,mode="markers"),surface_wavelet(wasg,200)]))
#display(PlotlyJS.plot([PlotlyJS.scatter(wasg,false,mode="markers"),PlotlyJS.surface(wasg,200)]))
#display(PlotlyJS.plot([PlotlyJS.scatter(wasg,false,mode="markers"),surface_wavelet(wasg,200)]))



import PlotlyJS
plots = PlotlyJS.GenericTrace{Dict{Symbol,Any}}[]
xx = linspace(-1,1,100)
for hcpt in wasg

	#if abs(scaling_weight(hcpt)) > .2#!isrefined(hcpt)
		#poly = x->scaling_weight(hcpt)*basis_fun(hcpt,[x],Maxp)
		#poly = x->hcpt.scaling_weight*basis_fun(hcpt,[x],Maxp)
		# poly = x->wavelet(wasg,hcpt,x],Maxp)
		display(PlotlyJS.plot(PlotlyJS.surface(poly, 20)))
		#println(getkey(hcpt)," => ",scaling_weight(hcpt))
		#push!(plots, PlotlyJS.scatter(x=xx, y=map(poly, xx)+level(hcpt), name="$(getkey(hcpt))") )
	#end
end
#PlotlyJS.plot(plots)



#for hcpt in wasg
#	if level(hcpt) > 1
#		println(hcpt)
#		fun = x->wavelet(hcpt,x,Maxp,Maxp-1)
#		display(PlotlyJS.plot(PlotlyJS.surface(fun,name="$(getkey(hcpt))")))
#	end
#end



##f =  x-> sqrt(sum(x.^2)) > .5 ? 1.0 : -1.0
##f =  x-> x[1] * x[2]> 2*pi/10. ? 1.0 : -1.0
#f =  x-> 1/(abs(sqrt(2)-(x[1]-1)^2-(x[2]-1)^2)+.5)
##f = x->  x[2]^2
#
#

#
#PlotlyJS.plot([PlotlyJS.scatter3d(wasg,mode="markers"),PlotlyJS.surface(wasg,400)])
##PlotlyJS.plot([PlotlyJS.scatter(wasg,false,mode="markers"),PlotlyJS.surface(wasg,1000)])

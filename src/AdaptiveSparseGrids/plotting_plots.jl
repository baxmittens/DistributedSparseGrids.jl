
using Plots
using Colors
import Colors: distinguishable_colors, RGB, N0f8, colormap

function scatterplot(sg::SG, lvl_offset::Float64=0.0, color_order::Bool=false, maxp::Int=1) where {CT,CP<:AbstractCollocationPoint{2,CT},HCP<:AbstractHierarchicalCollocationPoint{2,CP},SG<:AbstractHierarchicalSparseGrid{2,HCP}}
	if color_order
		colors = cols = colormap("Reds", N*maxp+1)
	else
		colors = cols = distinguishable_colors(numlevels(sg)+1, [RGB(1,1,1)])[2:end]
		#colors = cols = reverse(colormap("RdBu", numlevels(sg)+1))
		#colors = cols =  range(HSV(0,1,1), stop=HSV(-360,1,1), length=numlevels(sg)+1)
	end	

	nlevel = numlevels(sg)
	xvals = Vector{Vector{CT}}(undef,nlevel)
	yvals = Vector{Vector{CT}}(undef,nlevel)
	zvals = Vector{Vector{CT}}(undef,nlevel)
	text = Vector{Vector{String}}(undef,nlevel)
	clr = Vector{Vector{RGB{N0f8}}}(undef,nlevel)
	for l = 1:nlevel
		xvals[l] = Vector{CT}()
		yvals[l] = Vector{CT}()
		zvals[l] = Vector{CT}()
		clr[l] = Vector{RGB{N0f8}}()
		text[l] = Vector{String}()
	end
	for hcpt in sg
		l = level(hcpt)
		
		push!(xvals[l],coord(hcpt,1))
		push!(yvals[l],coord(hcpt,2))
		push!(zvals[l],level(hcpt)*lvl_offset)
		#push!(zvals,interpolate(sg, [xvals[end], yvals[end]]))
		push!(text[l],string(pt_idx(hcpt))*"^"*string(i_multi(hcpt)))
		if !color_order
			push!(clr[l],colors[level(hcpt)])
		else

			N = length(children(hcpt))
			ord = 0
			for d = 1:N
				ord += polyorder_v1(hcpt, d, maxp)
			end
			push!(clr[l],colors[ord+1])
		end
	end
	if lvl_offset > 0.0
		#p = PlotlyJS.scatter3d(x=xvals, y=yvals, z=zvals, text=text, marker_color=clr, mode="markers+text",marker_size=7,textposition="bottom center")
		#p = PlotlyJS.scatter3d(x=xvals[i], y=yvals[i], z=zvals[i], text=text[i], marker_color=clr[i], mode="markers",marker_size=2,textposition="bottom center")
		p = Plots.scatter(xvals[1], yvals[1], zvals[1])#, text=text[i], marker_color=clr[i], mode="markers",marker_size=mw,textposition="bottom center",name="level $i")
	else
		#p = PlotlyJS.scatter(x=xvals, y=yvals, text=text, marker_color=clr, mode="markers+text",marker_size=7,textposition="bottom center")
		#p = PlotlyJS.scatter(x=xvals[i], y=yvals[i], text=text[i], marker_color=clr[i], mode="markers",marker_size=2,textposition="bottom center")
		p = Plots.scatter(xvals[2], yvals[1])#, text=text[i], marker_color=clr[i], mode="markers",marker_size=mw,textposition="bottom center",name="level $i")
	end
	for i = 2:nlevel
		mw = 8.0-foldl((x,y)->x+2.0/(y),1:i)
		if lvl_offset > 0.0
			#p = PlotlyJS.scatter3d(x=xvals, y=yvals, z=zvals, text=text, marker_color=clr, mode="markers+text",marker_size=7,textposition="bottom center")
			#p = PlotlyJS.scatter3d(x=xvals[i], y=yvals[i], z=zvals[i], text=text[i], marker_color=clr[i], mode="markers",marker_size=2,textposition="bottom center")
			Plots.scatter!(p,xvals[i], yvals[i], zvals[i], color=clr[i])#, text=text[i], marker_color=clr[i], mode="markers",marker_size=mw,textposition="bottom center",name="level $i")
		else
			#p = PlotlyJS.scatter(x=xvals, y=yvals, text=text, marker_color=clr, mode="markers+text",marker_size=7,textposition="bottom center")
			#p = PlotlyJS.scatter(x=xvals[i], y=yvals[i], text=text[i], marker_color=clr[i], mode="markers",marker_size=2,textposition="bottom center")
			Plots.scatter!(p,xvals[i], yvals[i], color=clr[i])#, text=text[i], marker_color=clr[i], mode="markers",marker_size=mw,textposition="bottom center",name="level $i")
		end
	end
	return plot(p)
end

function scatterplot(sg::SG, lvl_offset::Bool=false; kwargs...) where {CT,CP<:AbstractCollocationPoint{1,CT},HCP<:AbstractHierarchicalCollocationPoint{1,CP},SG<:AbstractHierarchicalSparseGrid{1,HCP}}
	colors = cols = distinguishable_colors(numlevels(sg)+1, [RGB(1,1,1)])[2:end]
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	text = Vector{String}()
	clr = Vector{RGB{N0f8}}()
	for hcpt in sg
		push!(xvals,coord(hcpt,1))
		push!(yvals,level(hcpt))
		push!(text,string(pt_idx(hcpt))*"^"*string(i_multi(hcpt)))
		push!(clr,colors[level(hcpt)])
	end
	#p = Plots.scatter(xvals, lvl_offset ? yvals : zeros(CT,length(xvals)), text=text, color=clr)#; kwargs... )
	p = Plots.scatter(xvals, lvl_offset ? yvals : zeros(CT,length(xvals)), color=clr)#; kwargs... )
	return Plots.plot(p)
end

function surfaceplot(asg::SG, npts = 20, stoplevel::Int=numlevels(asg)) where {CT,CP<:AbstractCollocationPoint{2,CT},HCP<:AbstractHierarchicalCollocationPoint{2,CP},SG<:AbstractHierarchicalSparseGrid{2,HCP}}
	pts = range(-1.,stop=1.,length=npts)
	f(x,y) = interpolate(asg, [x,y], stoplevel)
	p = Plots.surface(pts,pts,f)
	return p
end

function lineplot(asg::SG, npts = 1000, stoplevel::Int=numlevels(asg); kwargs...) where {CT,CP<:AbstractCollocationPoint{1,CT},HCP<:AbstractHierarchicalCollocationPoint{1,CP},SG<:AbstractHierarchicalSparseGrid{1,HCP}}
	xpts = collect(range(-1.,stop=1.,length=npts))
	ypts = similar(xpts)
	for i = 1:npts
		ypts[i] = interpolate(asg, [xpts[i]], stoplevel)
	end
	p = Plots.plot(xpts,ypts; kwargs...)
	return p
end

function surfaceplot(fun::F, npts = 20; kwargs...) where {F<:Function}
	pts = range(-1.,stop=1.,length=npts)
	p = Plots.surface(pts,pts,fun)
	return p
end

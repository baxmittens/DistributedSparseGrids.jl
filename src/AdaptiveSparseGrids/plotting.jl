
using Colors
import Colors: distinguishable_colors, RGB, N0f8, colormap
import PlotlyJS
import PlotlyJS: GenericTrace
import UnicodePlots

function PlotlyJS.scatter(sg::SG, lvl_offset::Float64=0.0, color_order::Bool=false, maxp::Int=1) where {CT,CP<:AbstractCollocationPoint{2,CT},HCP<:AbstractHierarchicalCollocationPoint{2,CP},SG<:AbstractHierarchicalSparseGrid{2,HCP}}
	if color_order
		colors = cols = colormap("Reds", N*maxp+1)
	else
		colors = cols = distinguishable_colors(numlevels(sg)+1, [RGB(1,1,1)])[2:end]
		#colors = cols = reverse(colormap("RdBu", numlevels(sg)+1))
		#colors = cols =  range(HSV(0,1,1), stop=HSV(-360,1,1), length=numlevels(sg)+1)
	end	

	nlevel = numlevels(sg)
	traces = Vector{GenericTrace}(undef,nlevel)
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
	for i = 1:nlevel
		mw = 8.0-foldl((x,y)->x+2.0/(y),1:i)
		if lvl_offset > 0.0
			#p = PlotlyJS.scatter3d(x=xvals, y=yvals, z=zvals, text=text, marker_color=clr, mode="markers+text",marker_size=7,textposition="bottom center")
			#p = PlotlyJS.scatter3d(x=xvals[i], y=yvals[i], z=zvals[i], text=text[i], marker_color=clr[i], mode="markers",marker_size=2,textposition="bottom center")
			traces[i] = PlotlyJS.scatter3d(x=xvals[i], y=yvals[i], z=zvals[i], text=text[i], marker_color=clr[i], mode="markers",marker_size=mw,textposition="bottom center",name="level $i")
		else
			#p = PlotlyJS.scatter(x=xvals, y=yvals, text=text, marker_color=clr, mode="markers+text",marker_size=7,textposition="bottom center")
			#p = PlotlyJS.scatter(x=xvals[i], y=yvals[i], text=text[i], marker_color=clr[i], mode="markers",marker_size=2,textposition="bottom center")
			traces[i] = p = PlotlyJS.scatter(x=xvals[i], y=yvals[i], text=text[i], marker_color=clr[i], mode="markers",marker_size=mw,textposition="bottom center",name="level $i")
		end
	end
	return traces
end

function PlotlyJS.scatter(sg::SG, lvl_offset::Bool=false; kwargs...) where {CT,CP<:AbstractCollocationPoint{1,CT},HCP<:AbstractHierarchicalCollocationPoint{1,CP},SG<:AbstractHierarchicalSparseGrid{1,HCP}}
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
	p = PlotlyJS.scatter(x=xvals, y=lvl_offset ? yvals : zeros(CT,length(xvals)), text=text, marker_color=clr; kwargs... )
	return p
end

function UnicodePlots.scatterplot(sg::SG, lvl_offset::Bool=false) where {CT,CP<:AbstractCollocationPoint{1,CT},HCP<:AbstractHierarchicalCollocationPoint{1,CP},SG<:AbstractHierarchicalSparseGrid{1,HCP}}
	#colors = cols = distinguishable_colors(numlevels(sg)+1, [RGB(1,1,1)])[2:end]
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	#text = Vector{String}()
	#clr = Vector{RGB{N0f8}}()
	for hcpt in sg
		push!(xvals,coord(hcpt,1))
		push!(yvals,level(hcpt))
		#push!(text,string(pt_idx(hcpt))*"^"*string(i_multi(hcpt)))
		#push!(clr,colors[level(hcpt)])
	end
	#p = PlotlyJS.scatter(x=xvals, y=lvl_offset ? yvals : zeros(CT,length(xvals)), text=text, marker_color=clr; kwargs... )
	p = UnicodePlots.scatterplot(xvals, lvl_offset ? yvals : zeros(CT,length(xvals)))
	return p
end

function UnicodePlots.scatterplot(sg::SG) where {CT,CP<:AbstractCollocationPoint{2,CT},HCP<:AbstractHierarchicalCollocationPoint{2,CP},SG<:AbstractHierarchicalSparseGrid{2,HCP}}
	#colors = cols = distinguishable_colors(numlevels(sg)+1, [RGB(1,1,1)])[2:end]
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	#text = Vector{String}()
	#clr = Vector{RGB{N0f8}}()
	for hcpt in sg
		push!(xvals,coord(hcpt,1))
		push!(yvals,coord(hcpt,2))
		#push!(text,string(pt_idx(hcpt))*"^"*string(i_multi(hcpt)))
		#push!(clr,colors[level(hcpt)])
	end
	#p = PlotlyJS.scatter(x=xvals, y=lvl_offset ? yvals : zeros(CT,length(xvals)), text=text, marker_color=clr; kwargs... )
	p =  UnicodePlots.scatterplot(xvals, yvals)
	return p
end

function UnicodePlots.lineplot(asg::SG, npts = 1000, stoplevel::Int=numlevels(asg)) where {CT,CP<:AbstractCollocationPoint{1,CT},HCP<:AbstractHierarchicalCollocationPoint{1,CP},SG<:AbstractHierarchicalSparseGrid{1,HCP}}
	xpts = collect(range(-1.,stop=1.,length=npts))
	ypts = similar(xpts)
	for i = 1:npts
		yval = interpolate(asg, [xpts[i]], stoplevel)
		if typeof(yval) <: Number
			ypts[i] = yval
		else
			ypts[i] = norm(yval)
		end

	end
	p = UnicodePlots.lineplot(xpts,ypts)
end

function PlotlyJS.scatter3d(sg::SG, color_order::Bool=false, maxp::Int=1) where {CT,CP<:AbstractCollocationPoint{2,CT},HCP<:AbstractHierarchicalCollocationPoint{2,CP},SG<:AbstractHierarchicalSparseGrid{2,HCP}}
	if color_order
		colors = cols = colormap("Reds", N*maxp+1)
	else
		colors = cols = distinguishable_colors(numlevels(sg)+1, [RGB(1,1,1)])[2:end]
	end
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	zvals = Vector{CT}()
	text = Vector{String}()
	clr = Vector{RGB{N0f8}}()
	for hcpt in sg
		push!(xvals,coord(hcpt,1))
		push!(yvals,coord(hcpt,2))
		push!(zvals,zero(CT))
		push!(text,string(pt_idx(hcpt))*"^"*string(i_multi(hcpt)))
		if !color_order
			push!(clr,colors[level(hcpt)])
		else

			N = length(children(hcpt))
			ord = 0
			for d = 1:N
				ord += polyorder_v1(hcpt, d, maxp)
			end
			push!(clr,colors[ord+1])
		end
	end
	p = PlotlyJS.scatter3d(x=xvals, y=yvals, z=zvals, text=text, marker_color=clr, mode="markers",marker_size=2,textposition="bottom center")
	return p
end

function PlotlyJS.scatter3d(sg::SG, color_order::Bool=false, maxp::Int=1) where {CT,CP<:AbstractCollocationPoint{3,CT},HCP<:AbstractHierarchicalCollocationPoint{3,CP},SG<:AbstractHierarchicalSparseGrid{3,HCP}}
	if color_order
		colors = cols = colormap("Reds", N*maxp+1)
	else
		colors = cols = distinguishable_colors(numlevels(sg)+1, [RGB(1,1,1)])[2:end]
	end
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	zvals = Vector{CT}()
	text = Vector{String}()
	clr = Vector{RGB{N0f8}}()
	for hcpt in sg
		push!(xvals,coord(hcpt,1))
		push!(yvals,coord(hcpt,2))
		push!(zvals,coord(hcpt,3))
		push!(text,string(pt_idx(hcpt))*"^"*string(i_multi(hcpt)))
		if !color_order
			push!(clr,colors[level(hcpt)])
		else

			N = length(children(hcpt))
			ord = 0
			for d = 1:N
				ord += polyorder_v1(hcpt, d, maxp)
			end
			push!(clr,colors[ord+1])
		end
	end
	p = PlotlyJS.scatter3d(x=xvals, y=yvals, z=zvals, text=text, marker_color=clr, mode="markers",marker_size=2,textposition="bottom center")
	return p
end

include(joinpath("../support","ndgrid.jl"))
function PlotlyJS.surface(asg::SG, npts = 20, postfun=x->x; kwargs...) where {CT,CP<:AbstractCollocationPoint{2,CT},HCP<:AbstractHierarchicalCollocationPoint{2,CP},SG<:AbstractHierarchicalSparseGrid{2,HCP}}
	pts = range(-1.,stop=1.,length=npts)
	xpts, ypts = ndgrid(pts,pts)
	zz = similar(xpts)
	tmp = similar(scaling_weight(first(asg)))
	for i = 1:npts, j=1:npts
		fill!(tmp,0.0)
		interpolate!(tmp,asg, [xpts[i,j], ypts[i,j]])
		zz[i,j] = postfun(tmp)
	end
	p = PlotlyJS.surface(x=xpts,y=ypts,z=zz; kwargs...)
end

function surface_wavelet(asg::SG, npts = 20; kwargs...) where {CT,CP<:AbstractCollocationPoint{2,CT},HCP<:AbstractHierarchicalCollocationPoint{2,CP},SG<:AbstractHierarchicalSparseGrid{2,HCP}}
	pts = range(-1.,stop=1.,length=npts)
	xpts, ypts = ndgrid(pts,pts)
	zz = similar(xpts)
	for i = 1:npts, j=1:npts
		zz[i,j] = interpolate_wavelet(asg, [xpts[i,j], ypts[i,j]])
	end
	p = PlotlyJS.surface(x=xpts,y=ypts,z=zz; kwargs...)
end

function PlotlyJS.surface(asg::SG, npts = 1000, postfun=x->x, stoplevel::Int=numlevels(asg); kwargs...) where {CT,CP<:AbstractCollocationPoint{1,CT},HCP<:AbstractHierarchicalCollocationPoint{1,CP},SG<:AbstractHierarchicalSparseGrid{1,HCP}}
	xpts = collect(range(-1.,stop=1.,length=npts))
	ypts = similar(xpts)
	for i = 1:npts
		ypts[i] = postfun(interpolate(asg, [xpts[i]], stoplevel))
	end
	p = PlotlyJS.scatter(x=xpts,y=ypts; kwargs...)
end

function surface_wavelet(asg::SG, npts = 1000, stoplevel::Int=numlevels(asg) ; kwargs...) where {CT,CP<:AbstractCollocationPoint{1,CT},HCP<:AbstractHierarchicalCollocationPoint{1,CP},SG<:AbstractHierarchicalSparseGrid{1,HCP}}
	xpts = collect(range(-1.,stop=1.,length=npts))
	ypts = similar(xpts)
	for i = 1:npts
		ypts[i] = interpolate_wavelet(asg, [xpts[i]],stoplevel)
	end
	p = PlotlyJS.scatter(x=xpts,y=ypts; kwargs...)
end

function PlotlyJS.surface(fun::F, npts = 20; kwargs...) where {F<:Function}
	pts = range(-1.,stop=1.,length=npts)
	xpts, ypts = ndgrid(pts,pts)
	zz = similar(xpts)
	for i = 1:npts, j=1:npts
		zz[i,j] = fun([xpts[i,j], ypts[i,j]])
	end
	p = PlotlyJS.surface(x=xpts,y=ypts,z=zz; kwargs...)
end

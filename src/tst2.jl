include("AdaptiveSparseGrids.jl")
include("testfuncs.jl")

function test_wavelet(N,maxp,f,tol,int_exact,nrefsteps,maxlvl,nvmom=maxp)
	CT = Float64
	RT = Float64
	CPType = CollocationPoint{N,CT}

	HCPType = HierarchicalCollocationPoint{N,CPType,RT}

	Maxp = maxp
	Nv = numwaveletcoeffs(N, nvmom)
	force_zero_mean = true
	WCPType = WaveletCollocationPoint{N,CPType,RT,CT,Nv}

	pointprobs = SVector{N,Int}([1 for i = 1:N])

	wasg = init(AHSG{N,WCPType},pointprobs,Maxp)
	for i = 1:4; generate_next_level!(wasg); end
	wavelet_coeffs!(wasg, nvmom)
	init_weights!(wasg, f)
	wavelet_weights!(wasg,f,nvmom)

	errs = Vector{Float64}()
	#ints = Vector{Float64}()
	npts = Vector{Int}()
	push!( errs, abs(integrate(wasg)-int_exact)/int_exact)
	#push!( ints, integrate(wasg) )
	push!(npts, numpoints(wasg))

	for i = 1:nrefsteps
		@info "Refinement step $i/$nrefsteps"
		nchilds = wavelet_refine!(wasg,tol,maxlvl)
		@info "$(length(nchilds)) new CollocationPoints created"
		if !isempty(nchilds)
			wavelet_coeffs!(wasg, nvmom)
			init_weights!(wasg, f)
			wavelet_weights!(wasg,f,nvmom)
			#push!( errs, abs(abs(int_exact-integrate(wasg))/int_exact) )
			interr = abs(integrate(wasg)-int_exact)/int_exact
			println("$i npts = $(npts[end]), interr = $interr")
			push!( errs,  interr )
			push!(npts, numpoints(wasg))
		end
	end

	return wasg,PlotlyJS.scatter(x=npts,y=errs,name="wavlet_asg"),PlotlyJS.scatter(wasg,false)
end

function test_lagrange(N,maxp,f,tol,int_exact,nrefsteps,maxlvl)
	CT = Float64
	RT = Float64
	CPType = CollocationPoint{N,CT}
	HCPType = HierarchicalCollocationPoint{N,CPType,RT}
	Maxp = maxp
	pointprobs = SVector{N,Int}([1 for i = 1:N])
	wasg = init(AHSG{N,HCPType},pointprobs,Maxp)
	for i = 1:6; generate_next_level!(wasg); end
	init_weights!(wasg, f)

	errs = Vector{Float64}()

	npts = Vector{Int}()

	#push!( errs, abs(integrate(wasg)-pi)/pi)
	#push!(npts, numpoints(wasg))

	for i = 1:nrefsteps
		@info "refstep $i"
		@time begin
			nchilds = generate_next_level!(wasg,tol,maxlvl)
			@info "$(length(nchilds)) new CollocationPoints created"
			if !isempty(nchilds)
	
				@info "start init weights"
				t1 = time()
				#init_weights!(wasg, f)
				init_weights!(wasg, nchilds, f) 
				t2 = time()
				println(t2-t1,"s init weights")
				
				
				#if numpoints(wasg)>300 && numpoints(wasg)<3000
					@info "store error"
					t1 = time()
					push!( errs, abs(integrate(wasg)-int_exact)/int_exact)
					#push!( errs, abs(integrate(wasg)))
					t2 = time()
					println(t2-t1,"s integrate")
					push!(npts, numpoints(wasg))
				#end
			end
		end
	end

	return wasg,PlotlyJS.scatter(x=npts,y=errs,name="Sparse Grid"),PlotlyJS.scatter(wasg,0.0,false)
end

import PlotlyJS

N=2
maxp=1
x0 = [pi/20 for i = 1:N]
#<<<<<<< HEAD
#f = x -> NBALL(x)
f =  x -> 1/(abs(sqrt(2)-(x[1]-1)^2-(x[2]-1)^2)+.5)
#@time mcint = MC_integral(Val{N},Val{1_000_000_000},f)
int_exact = BigFloat(2.69700031935962954830746595863370355419148701026617238036194035998592255712069)
nrefsteps = 20
maxlvl = 8
plots = PlotlyJS.GenericTrace{Dict{Symbol,Any}}[]
tol = 1e-16
#wasg,p1,ptsplot1 = test_wavelet(N,maxp,f,tol,int_exact,nrefsteps,maxlvl,maxp)
# =======
# #f = x -> NBALL(x)
# f =  x -> 1/(abs(sqrt(2)-(x[1]-1)^2-(x[2]-1)^2)+.5)
# @time mcint = MC_integral(Val{N},Val{100_000_000},f)
# int_exact = mcint#pi
# nrefsteps = 20
# plots = PlotlyJS.GenericTrace{Dict{Symbol,Any}}[]
# tol = .05
# wasg,p1,ptsplot1 = test_wavelet(N,maxp,f,tol,int_exact,nrefsteps,maxp-1)
# >>>>>>> df4103520d995b4be09632feb2141fa244210cf3
#push!(plots,p1)
_hasg,_p2,_ptsplot2 = test_lagrange(N,maxp,f,tol,int_exact,nrefsteps,maxlvl)
if 0==1
@info "$(numpoints(hasg)) collocation point"
push!(plots,p2)

import PlotlyJS: attr, Layout
 axis = attr(showgrid=false, gridcolor="white", linewidth=1.0,
                linecolor="white", titlefont_color="#555555",
                titlefont_size=14, showticks=false,
                tickcolor="#555555",showticklabels=false,zeroline=false,showline=false,tickwidth=0,ticklen=0
                )
 #ax = attr( backgroundcolor="#FF5555",gridcolor="white",showbackground=true,nticks=5,tickfont_size=12)
 ax = attr( nticks=5,tickfont_size=12)
 zax = attr( nticks=5,tickfont_size=12,title="level")
 ax1 = attr(scaleanchor="x",scalaratio = 1.0)
 #PlotlyJS.plot(PlotlyJS.scatter(hasg,0.25),PlotlyJS.Layout(plot_bgcolor="white", paper_bgcolor="white",xaxis=axis,yaxis=axis,zaxis=axis))
 #p=PlotlyJS.plot(PlotlyJS.scatter(hasg),Layout(yaxis=ax1))
# p=PlotlyJS.plot(PlotlyJS.scatter(hasg,1.00),PlotlyJS.Layout(font_size=20,scene=attr(xaxis=ax,yaxis=ax,zaxis=zax,aspectmode="manual",aspectratio=attr(x=2, y=2, z=1))))
#PlotlyJS.plot(plots,PlotlyJS.Layout(xaxis_type="log",yaxis_type="log"))
using Statistics

res = Vector{Vector{Float64}}()
vars1 = Float64[]
vars2 = Float64[]
boxes = GenericTrace[]
numpts = [10000,100000,1000000]#,1000000]#,10000000,100000000]
for numpt in numpts
	@info "$numpt points"
	intres = Float64[]
	#for i = 1:100
	for i = 1:1000
		push!(intres,MC_integral(Val{N},Val{numpt},f))
		#if numpt >= 1000000
		#	break
		#end
	end
	push!(res,intres)
end

#xx = vcat(numpts,reverse(numpts))
#yy = [abs(mean(res[1])-int_exact+var(res[1]))/int_exact,abs(mean(res[2])-int_exact+var(res[2]))/int_exact,abs(mean(res[3])-int_exact+var(res[3]))/int_exact,abs(mean(res[1])-int_exact-var(res[1]))/int_exact,abs(mean(res[2])-int_exact-var(res[2]))/int_exact,abs(mean(res[3])-int_exact-var(res[3]))/int_exact]

#trace1 = scatter(x=xx,y=yy,fill="tozerox",fillcolor="rgba(0, 100, 80, 0.2)",line_color="transparent",name="Fair",showlegend=false)
#trace2 = scatter(x=numpts,y=[abs(mean(res[1])-int_exact)/int_exact,abs(mean(res[2])-int_exact)/int_exact,abs(mean(res[3])-int_exact)/int_exact],name="Monte Carlo") 
#trace3 = box(x=[numpts[1] for i = 1:length(res[1])],y=abs.(res[1].-int_exact)./int_exact,boxmean=true)
#trace4 = box(x=[numpts[2] for i = 1:length(res[2])],y=abs.(res[2].-int_exact)./int_exact,boxmean=true)
#trace5 = box(x=[numpts[3] for i = 1:length(res[3])],y=abs.(res[3].-int_exact)./int_exact,boxmean=true)

#push!(plots,trace1)
#traceMC = scatter(x=numpts,y=map(x->abs(mean(x.-int_exact))/int_exact,res),name="Monte Carlo")
#tracesMC = map((y,x)->PlotlyJS.scatter(x=[x for i = 1:length(y)], y=y, mode="markers") ,res, numpts)
#tracesMC = map((y,x)->PlotlyJS.box(x=[x for i = 1:length(y)], y=abs.(y.-int_exact)/int_exact, boxpoints=false, boxmean=false, name="Monte Carlo",marker_size=1,jitter=0.1,whiskerwidth=0.1) ,res, numpts)
#nmp = vcat(map((y,x)->[x for i = 1:length(y)],res,numpts)...)
#vcres = (vcat(res...).-int_exact)./int_exact
#tracesMC=PlotlyJS.box(x=nmp, y=vcres, boxpoints=false, boxmean=false, name="Monte Carlo",marker_size=1,whiskerwidth=0.2)
fz=18
xx = vcat(map((x,y)->[y for i = 1:length(x)],res,numpts)...)
tracesMC=PlotlyJS.box(x=xx,y=abs.(vcat(res...).-int_exact)./int_exact,mode="markers",boxpoints=false,name="Monte Carlo")
layout = Layout(xaxis_type = "log",yaxis_type = "log",yaxis=attr(automargin=true,title=attr(text="%error",standoff=2)),xaxis=attr(automargin=true,title=attr(text="#points",standoff=1)),showlegend=true,legend=attr(x=0.325,y=1.05,orientation="h",achor="center"),font_size=fz)
p=PlotlyJS.plot(vcat(plots,tracesMC),layout)#,height=500,width=800,font_size=fz))
end
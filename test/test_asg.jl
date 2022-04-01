include(joinpath("/home","bittens","workspace","AdaptiveSparseGrids","src","AdaptiveSparseGrids.jl"))
include(joinpath("/home","bittens","workspace","AdaptiveSparseGrids","src","testfuncs.jl"))

using UnicodePlots
fun(x,ID) = testfunc(x)
N=5
CT = Float64
RT = Float64
CPType = CollocationPoint{N,CT}
HCPType = HierarchicalCollocationPoint{N,CPType,RT}

Maxp = 1
maxlvl = 20
nrefsteps = 10
tol = 1e-5
pointprobs = SVector{N,Int}([1 for i = 1:N])
wasg = init(AHSG{N,HCPType},pointprobs,Maxp)
_cpts = Set{HierarchicalCollocationPoint{N,CPType,RT}}(collect(wasg))
for i = 1:2; union!(_cpts,generate_next_level!(wasg)); end
init_weights!(wasg, _cpts, fun)

for i = 1:nrefsteps
	@info "refstep $i"
	@time nchilds = generate_next_level!(wasg,tol,maxlvl)
	@info "$(length(nchilds)) new CollocationPoints created"
	if !isempty(nchilds)
		@info "start init weights"
		@time init_weights!(wasg, nchilds, fun)  
	end
end
using DistributedSparseGrids
import DistributedSparseGrids: refine!
using StaticArrays 


function getnextcpts(asg)
	nl = DistributedSparseGrids.numlevels(asg)
	cpts = filter(x->DistributedSparseGrids.level(x)==nl,collect(asg))
	if nl > 1
		filter!(cpt->all(cpt.cpt.coords .== 0.0 .|| cpt.cpt.coords .== 1.0 .|| cpt.cpt.coords .== -1.0), cpts)
	end
	return cpts
end

function refineedges!(asg, nlvl=3)
	for i = 1:nlvl-1
		cpts = getnextcpts(asg)
		map(x->refine!(asg,x),cpts)
	end
end

function gethyperedges(asg::DistributedSparseGrids.AdaptiveHierarchicalSparseGrid{N}) where N
	cpts = filter(x->DistributedSparseGrids.level(x)==N+1,collect(asg))
	if nl > 1
		filter!(cpt->all(cpt.cpt.coords .== 1.0 .|| cpt.cpt.coords .== -1.0), cpts)
	end
	return cpts
end



N = 9
RT=Float64
CT=Float64
CPType = CollocationPoint{N,CT}
HCPType = HierarchicalCollocationPoint{N,CPType,RT}
pointprops = @SVector ones(Int,N)
asg = init(AHSG{N,HCPType},pointprops)
refineedges!(asg, 10)
gethyperedges(asg)
DistributedSparseGrids.numpoints(asg)

cpts = getnextcpts(asg)
map(x->refine!(asg,x),cpts)





cpts = collect(asg)
lvl2cpts = refine!(asg, cpts[1])
cpts = filter(x->DistributedSparseGrids.level(x)==2,collect(asg))
lvl3cpts = refine!(asg, cpts[1])
cpts = filter(x->DistributedSparseGrids.level(x)==3,collect(asg))
lvl4cpts = refine!(asg, cpts[1])

cpts = filter(x->DistributedSparseGrids.level(x)==4,collect(asg))
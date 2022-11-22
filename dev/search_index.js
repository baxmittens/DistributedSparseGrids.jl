var documenterSearchIndex = {"docs":
[{"location":"lib/lib/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"lib/lib/#Contents","page":"Library","title":"Contents","text":"","category":"section"},{"location":"lib/lib/","page":"Library","title":"Library","text":"Pages = [\"lib.md\"]\nDepth = 4","category":"page"},{"location":"lib/lib/#Functions","page":"Library","title":"Functions","text":"","category":"section"},{"location":"lib/lib/#Index","page":"Library","title":"Index","text":"","category":"section"},{"location":"lib/lib/","page":"Library","title":"Library","text":"Pages = [\"lib.md\"]","category":"page"},{"location":"lib/lib/#Typedefs","page":"Library","title":"Typedefs","text":"","category":"section"},{"location":"lib/lib/","page":"Library","title":"Library","text":"DistributedSparseGrids.AbstractSparseGrid\nDistributedSparseGrids.AbstractHierarchicalSparseGrid\nDistributedSparseGrids.PointDict","category":"page"},{"location":"lib/lib/#DistributedSparseGrids.AbstractSparseGrid","page":"Library","title":"DistributedSparseGrids.AbstractSparseGrid","text":"AbstractSparseGrid{N}\n\nAbstract Type\n\n`N` : Dimension of hierarchical sparse grid\n\n\n\n\n\n","category":"type"},{"location":"lib/lib/#DistributedSparseGrids.AbstractHierarchicalSparseGrid","page":"Library","title":"DistributedSparseGrids.AbstractHierarchicalSparseGrid","text":"AbstractHierarchicalSparseGrid{N,HCP}\n\nAbstract Type\n\n`N` : Dimension of hierarchical sparse grid\n`HCP<:AbstractHierarchicalCollocationPoint` : Collocation point type\n\n\n\n\n\n","category":"type"},{"location":"lib/lib/#DistributedSparseGrids.PointDict","page":"Library","title":"DistributedSparseGrids.PointDict","text":"PointDict{ N, HCP <: AbstractHierarchicalCollocationPoint{N}}\n\nTypedef for `Dict{SVector{N,Int},Dict{SVector{N,Int},HCP}}`\n\n`N` : Dimension of hierarchical sparse grid\n`HCP<:AbstractHierarchicalCollocationPoint` : Collocation point type\n\n\n\n\n\n","category":"type"},{"location":"lib/lib/#Structs","page":"Library","title":"Structs","text":"","category":"section"},{"location":"lib/lib/","page":"Library","title":"Library","text":"DistributedSparseGrids.AdaptiveHierarchicalSparseGrid","category":"page"},{"location":"lib/lib/#DistributedSparseGrids.AdaptiveHierarchicalSparseGrid","page":"Library","title":"DistributedSparseGrids.AdaptiveHierarchicalSparseGrid","text":"AdaptiveHierarchicalSparseGrid{N,HCP}\n\nContainer for hierarchical collocation points \n\n# Fields\n\n`cpts::Vector{PointDict{N,HCP}}` : [`DistributedSparseGrids.PointDict`](@ref) with collocation points\n`pointSetProperties::SVector{N,Int}` : Vector containing all pointset properties. \n\nPointset properties = [psp_1,...,psp_N], \npsp_i in [1,2,3,4]. \n1=>`closed point set`, \n2=>`open point set`, \n3=>`left-open point set`, \n4=>`right-open point set`.\n\n\n\n\n\n","category":"type"},{"location":"lib/lib/#General-functions","page":"Library","title":"General functions","text":"","category":"section"},{"location":"lib/lib/","page":"Library","title":"Library","text":"init\ninterpolate\nintegrate\ninit_weights!\ninit_weights_inplace_ops!\ndistributed_init_weights!\ndistributed_init_weights_inplace_ops!","category":"page"},{"location":"lib/lib/#DistributedSparseGrids.init","page":"Library","title":"DistributedSparseGrids.init","text":"init(::Type{AHSG{N,HCP}}, pointSetProperties::SVector{N,Int})\n\nInitialize the sparse grid. Returns a N-dimensional sparse grid where only the root point has been created.\n\nConstructor\n\n::Type{AHSG{N,HCP}}: Define type of DistributedSparseGrids.ASHG\npointSetProperties::SVector{N,Int}: Vector containing all pointset properties.\n\nPointset properties = [psp1,...,pspN], psp_i in [1,2,3,4].  1=>closed point set,  2=>open point set,  3=>left-open point set,  4=>right-open point set.\n\nExample\n\nN = 1 pointprobs = @SVector [1] RT = Float64 CT = Float64 CPType = CollocationPoint{N,CT} HCPType = HierarchicalCollocationPoint{N,CPType,RT} asg = init(AHSG{N,HCPType},pointprobs)\n\n\n\n\n\n","category":"function"},{"location":"lib/lib/#DistributedSparseGrids.interpolate","page":"Library","title":"DistributedSparseGrids.interpolate","text":"interpolate(asg::SG, x::VCT, stplvl::Int=numlevels(asg))\n\nInterpolate at position x. \n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: initialized adaptive sparse grid\ncpts::Dict{SVector{N,Int},Dict{SVector{N,Int},HCP}}: Dict cointaining all collocation points\n\n\n\n\n\n","category":"function"},{"location":"lib/lib/#DistributedSparseGrids.integrate","page":"Library","title":"DistributedSparseGrids.integrate","text":"integrate(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}\n\nComputes the VTUHeader based on the headertype and a Base64 decoded input data array.\n\nConstructor\n\n::Type{T}: headertype, either UInt32 or UInt64\ninput::Vector{UInt8}: input data\n\nFields\n\nnum_blocks::T : number of blocks\nblocksize::T : size of blocks\nlast_blocksize::T : size of last block (can be different)\ncompressed_blocksizes::T : size of compressed blocks\n\n\n\n\n\n","category":"function"},{"location":"lib/lib/#DistributedSparseGrids.init_weights!","page":"Library","title":"DistributedSparseGrids.init_weights!","text":"init_weights!(asg::SG, cpts::AbstractVector{HCP}, fun::F)\n\nComputes all weights in cpts. \n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: adaptive sparse grid\ncpts::AbstractVector{HCP}: all weights of the collocation points in cpts will be (re-)calculated.\nfun::Function to be interpolated.\n\n\n\n\n\ninit_weights!(asg::SG, fun::F)\n\nComputes all weights in asg. \n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: adaptive sparse grid\nfun::Function to be interpolated.\n\n\n\n\n\n","category":"function"},{"location":"lib/lib/#DistributedSparseGrids.init_weights_inplace_ops!","page":"Library","title":"DistributedSparseGrids.init_weights_inplace_ops!","text":"init_weights_inplace_ops!(asg::SG, cpts::AbstractVector{HCP}, fun::F)\n\n(Re-)Computes all weights in cpts with in-place operations. In-place functions mul!(::RT,::RT),mul!(::RT,::Float64),add!(::RT,::RT) have to be defined.\n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: adaptive sparse grid\ncpts::AbstractVector{HCP}: all weights of the collocation points in cpts will be (re-)calculated.\nfun::Function to be interpolated.\n\n\n\n\n\ninit_weights_inplace_ops!(asg::SG, fun::F)\n\n(Re-)Computes all weights in asg with in-place operations. In-place functions mul!(::RT,::RT),mul!(::RT,::Float64),add!(::RT,::RT) have to be defined.\n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: adaptive sparse grid\n\n\n\n\n\n","category":"function"},{"location":"lib/lib/#DistributedSparseGrids.distributed_init_weights!","page":"Library","title":"DistributedSparseGrids.distributed_init_weights!","text":"distributed_init_weights!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int})\n\nComputes all weights in cpts on all workers in worker_ids. \n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: adaptive sparse grid\ncpts::AbstractVector{HCP}: all weights of the collocation points in cpts will be (re-)calculated.\nfun::Function to be interpolated.\nworker_ids: All available workers (can be added via using Distributed; addprocs(...)).\n\n\n\n\n\ndistributed_init_weights!(asg::SG, fun::F, worker_ids::Vector{Int})\n\nComputes all weights in asg on all workers in worker_ids. \n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: adaptive sparse grid\nfun::Function to be interpolated.\nworker_ids: All available workers (can be added via using Distributed; addprocs(...)).\n\n\n\n\n\n","category":"function"},{"location":"lib/lib/#DistributedSparseGrids.distributed_init_weights_inplace_ops!","page":"Library","title":"DistributedSparseGrids.distributed_init_weights_inplace_ops!","text":"distributed_init_weights_inplace_ops!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int})\n\nComputes all weights in cpts on all workers in worker_ids. In-place functions mul!(::RT,::RT),mul!(::RT,::Float64),add!(::RT,::RT) have to be defined.\n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: adaptive sparse grid\ncpts::AbstractVector{HCP}: all weights of the collocation points in cpts will be (re-)calculated.\nfun::Function to be interpolated.\nworker_ids: All available workers (can be added via using Distributed; addprocs(...)).\n\n\n\n\n\ndistributed_init_weights_inplace_ops!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int})\n\nComputes all weights in cpts on all workers in worker_ids. In-place functions mul!(::RT,::RT),mul!(::RT,::Float64),add!(::RT,::RT) have to be defined.\n\nArguments\n\nasg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}: adaptive sparse grid\ncpts::AbstractVector{HCP}: all weights of the collocation points in cpts will be (re-)calculated.\nfun::Function to be interpolated.\nworker_ids: All available workers (can be added via using Distributed; addprocs(...)).\n\n\n\n\n\n","category":"function"},{"location":"lib/lib/#Utils","page":"Library","title":"Utils","text":"","category":"section"},{"location":"lib/lib/","page":"Library","title":"Library","text":"DistributedSparseGrids.AHSG","category":"page"},{"location":"lib/lib/#DistributedSparseGrids.AHSG","page":"Library","title":"DistributedSparseGrids.AHSG","text":"const AHSG{N,HCP} = AdaptiveHierarchicalSparseGrid{N,HCP}\n\nShortcut for AdaptiveHierarchicalSparseGrid\n\n\n\n\n\n","category":"type"},{"location":"#DistributedSparseGrids.jl","page":"Home","title":"DistributedSparseGrids.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia library that implements an Adaptive Sparse Grid collocation method for integrating memory-heavy objects generated on distributed workers (link to GitHub repository).","category":"page"},{"location":"","page":"Home","title":"Home","text":"For an alternative implementation, see AdaptiveSparseGrids.jl.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"index.md\", \"lib/lib.md\"]\nDepth = 3","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To mitigate the \"curse of dimensionality\" that occurs in the integration or interpolation of high-dimensional functions using tensor-product discretizations, sparse grids use Smolyak's quadrature rule. This is particularly useful if the evaluation of the underlying function is costly. In this library, an Adaptive Sparse Grid Collocation method with a local hierarchical Lagrangian basis, first proposed by Ma and Zabaras (2010), is implemented. For more information about the construction of Sparse Grids, see e.g. Gates and Bittens (2015).","category":"page"},{"location":"#Install","page":"Home","title":"Install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"import Pkg\nPkg.install(\"DistributedSparseGrids\")","category":"page"},{"location":"#Implemented-features","page":"Home","title":"Implemented features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"-\tNested one-dimensional Clenshaw-Curtis rule","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\tSmolyak's sparse grid construction","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\tlocal hierarchical Lagrangian basis","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\tdifferent pointsets (open, closed, halfopen)","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\tadaptive refinement","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\tdistributed function evaluation with Distributed.remotecall_fetch","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\tmulti-threaded calculation of basis coefficients with Threads.@threads","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\tusage of arbitrary return types ","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\tintegration","category":"page"},{"location":"","page":"Home","title":"Home","text":"-\texperimental: integration over X_sim (i) (the X_sim (i)  notation indicates the set of all variables except X_i).","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"#Point-sets","page":"Home","title":"Point sets","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"When using sparse grids, one can choose whether the 2d second-level collocation points should lay on the boundary of the domain or in the middle between the origin and the boundary. (There are other choices as well.) This results in two different sparse grids, the former with almost all points on the boundary and on the coordinate axes, the latter with all points in the interior of the domain. Since one can choose for both one-dimensional children of the root point individually, there exist a multitude of different point sets for Sparse Grids.","category":"page"},{"location":"","page":"Home","title":"Home","text":"DistributedSparseGrids\nusing StaticArrays \n\nfunction sparse_grid(N::Int,pointprobs,nlevel=6,RT=Float64,CT=Float64)\n\t# define collocation point\n\tCPType = CollocationPoint{N,CT}\n\t# define hierarchical collocation point\n\tHCPType = HierarchicalCollocationPoint{N,CPType,RT}\n\t# init grid\n\tasg = init(AHSG{N,HCPType},pointprobs)\n\t#set of all collocation points\n\tcpts = Set{HierarchicalCollocationPoint{N,CPType,RT}}(collect(asg))\n\t# fully refine grid nlevel-1 times\n\tfor i = 1:nlevel-1\n\t\tunion!(cpts,generate_next_level!(asg))\n\tend\n\treturn asg\nend\n\n# define point properties \n#\t1->closed point set\n# \t2->open point set\n#\t3->left-open point set\n#\t4->right-open point set\n\nasg01 = sparse_grid(1, @SVector [1]) \nasg02 = sparse_grid(1, @SVector [2]) \nasg03 = sparse_grid(1, @SVector [3]) \n\nasg04 = sparse_grid(2, @SVector [1,1]) \nasg05 = sparse_grid(2, @SVector [2,2]) \nasg06 = sparse_grid(2, @SVector [1,2]) \nasg07 = sparse_grid(2, @SVector [2,1]) \nasg08 = sparse_grid(2, @SVector [3,3]) \nasg09 = sparse_grid(2, @SVector [4,4]) \nasg10 = sparse_grid(2, @SVector [3,1]) \nasg11 = sparse_grid(2, @SVector [2,3]) \nasg12 = sparse_grid(2, @SVector [4,2]) ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: ) (Image: )","category":"page"},{"location":"#Integration-and-Interpolation","page":"Home","title":"Integration and Interpolation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"asg = sparse_grid(4, @SVector [1,1,1,1]) \n\n#define function: input are the coordinates x::SVector{N,CT} and an unique id ID::String (e.g. \"1_1_1_1\")\nfun1(x::SVector{N,CT},ID::String) = sum(x.^2)\n\n# initialize weights\n@time init_weights!(asg, fun1)\n\n# integration\nintegrate(asg)\n\n# interpolation\nx = rand(4)*2.0 .- 1.0\nval = interpolate(asg,x)\t","category":"page"},{"location":"#Distributed-function-evaluation","page":"Home","title":"Distributed function evaluation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"asg = sparse_grid(4, @SVector [1,1,1,1]) \n\n# add worker and register function to all workers\nusing Distributed\naddprocs(2)\nar_worker = workers()\n@everywhere begin\n    using StaticArrays\n    fun2(x::SVector{4,Float64},ID::String) = 1.0\nend\n\n# Evaluate the function on 2 workers\ndistributed_init_weights!(asg, fun2, ar_worker)","category":"page"},{"location":"#Using-custom-return-types","page":"Home","title":"Using custom return types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For custom return type T to work, following functions have to be implemented","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Base: +,-,*,/,^,zero,zeros,one,ones,copy,deepcopy\n\n+(a::T, b::T) \n+(a::T, b::Float64) \n*(a::T, b::Float64) \n-(a::T, b::Matrix{Float64})\n-(a::T, b::Float64) \nzero(a::T) \nzeros(a::T) \none(a::T) \none(a::T) \ncopy(a::T)\ndeepcopy(a::T)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is already the case for many data types. Below  RT=Matrix{Float64} is used.","category":"page"},{"location":"","page":"Home","title":"Home","text":"# sparse grid with 5 dimensions and levels\npointprop = @SVector [1,2,3,4,1]\nasg = sparse_grid(5, pointprop, 6, Matrix{Float64}) \n\n# define function: input are the coordinates x::SVector{N,CT} and an unique id ID::String (e.g. \"1_1_1_1_1_1_1_1_1_1\"\n# for the root poin in five dimensions)\nfun3(x::SVector{N,CT},ID::String) = ones(100,100).*x[1]\n\n# initialize weights\n@time init_weights!(asg, fun3)","category":"page"},{"location":"#In-place-operations","page":"Home","title":"In-place operations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are many mathematical operations executed which allocate memory while evaluting the hierarchical interpolator. Many of these allocations can be avoided by additionally implementing the inplace operations interface for data type T.","category":"page"},{"location":"","page":"Home","title":"Home","text":"import LinearAlgebra\nimport LinearAlgebra: mul!\n\nDistributedSparseGrids.add!(a::T, b::T) \nDistributedSparseGrids.add!(a::T, b::Float64) \nLinearAlgebra.mul!(a::T, b::Float64) \nLinearAlgebra.mul!(a:T, b::T, c::Float64)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For Matrix{Float64} this interface is already implemented.","category":"page"},{"location":"","page":"Home","title":"Home","text":"# initialize weights\n@time init_weights_inplace_ops!(asg, fun3)","category":"page"},{"location":"#Distributed-function-evaluation-and-in-place-operations","page":"Home","title":"Distributed function evaluation and in-place operations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"# initialize weights\n@time distributed_init_weights_inplace_ops!(asg, fun3, ar_worker)","category":"page"},{"location":"#Adaptive-Refinement","page":"Home","title":"Adaptive Refinement","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"# Sparse Grid with 4 initial levels\npp = @SVector [1,1]\nasg = sparse_grid(2, pp, 4)\n\n# Function with curved singularity\nfun1(x::SVector{2,Float64},ID::String) =  (1.0-exp(-1.0*(abs(2.0 - (x[1]-1.0)^2.0 - (x[2]-1.0)^2.0) +0.01)))/(abs(2-(x[1]-1.0)^2.0-(x[2]-1.0)^2.0)+0.01)\n\ninit_weights!(asg, fun1)\n\n# adaptive refine\nfor i = 1:20\n# call generate_next_level! with tol=1e-5 and maxlevels=20\ncpts = generate_next_level!(asg, 1e-5, 20)\ninit_weights!(asg, collect(cpts), fun1)\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Plotting","page":"Home","title":"Plotting","text":"","category":"section"},{"location":"#d","page":"Home","title":"1d","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"# grid plots\nPlotlyJS.scatter(sg::AbstractHierarchicalSparseGrid{1,HCP}, lvl_offset::Bool=false; kwargs...) \nUnicodePlots.scatterplot(sg::AbstractHierarchicalSparseGrid{1,HCP}, lvl_offset::Bool=false)\n\n# response function plots\nUnicodePlots.lineplot(asg::AbstractHierarchicalSparseGrid{1,HCP}, npts = 1000, stoplevel::Int=numlevels(asg))\nPlotlyJS.surface(asg::SG, npts = 1000, stoplevel::Int=numlevels(asg); kwargs...)","category":"page"},{"location":"#d-2","page":"Home","title":"2d","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"# grid plots\nPlotlyJS.scatter(sg::AbstractHierarchicalSparseGrid{2,HCP}, lvl_offset::Float64=0.0, color_order::Bool=false) \nUnicodePlots.scatterplot(sg::AbstractHierarchicalSparseGrid{2,HCP})\nPlotlyJS.scatter3d(sg::AbstractHierarchicalSparseGrid{2,HCP}, color_order::Bool=false, maxp::Int=1)\n\n# response function plot\nPlotlyJS.surface(asg::AbstractHierarchicalSparseGrid{2,HCP}, npts = 20; kwargs...)","category":"page"},{"location":"#d-3","page":"Home","title":"3d","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"# grid plot\nPlotlyJS.scatter3d(sg::AbstractHierarchicalSparseGrid{3,HCP}, color_order::Bool=false, maxp::Int=1)","category":"page"}]
}

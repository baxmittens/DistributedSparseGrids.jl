import Pkg; Pkg.add("Documenter")
push!(LOAD_PATH,"../src/")
using Documenter, DistributedSparseGrids
makedocs(
	sitename = "DistributedSparseGrids.jl",
	modules = [DistributedSparseGrids],
	pages = [
		"Home" => "index.md"
		"Library" => "lib/lib.md"
	]
	)
deploydocs(
    repo = "github.com/baxmittens/DistributedSparseGrids.jl.git"
)
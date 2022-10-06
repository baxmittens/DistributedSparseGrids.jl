# VTUFileHandler.jl
A VTU library in the Julia language that implements an algebra for basic mathematical operations on VTU data ([link to GitHub repository](https://github.com/baxmittens/VTUFileHandler.jl)).

## Contents

```@contents
Pages = ["index.md", "lib/lib.md"]
Depth = 2
```

## Introduction

With increasing computing resources, investigating uncertainties in simulation results is becoming an increasingly important factor. A discrete numerical simulation is computed several times with different deviations of the input parameters to produce different outputs of the same model to analyze those effects. The relevant stochastic or parametric output variables, such as mean, expected value, and variance, are often calculated and visualized only at selected individual points of the whole domain. This project aims to provide a simple way to perform stochastic/parametric post-processing of numerical simulations on entire domains using the [VTK unstructured grid](https://vtk.org/) (VTU) [file system](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) and the Julia language as an example.

## Install

```julia
import Pkg
Pkg.add("VTUFileHandler")
```

## Preliminaries 

The VTUFileHandler will eventually be used to perform stochastic post-processing on large VTU result files. Therefore, the following assumptions have to be fulfilled for the software to work correctly:

1. The VTU file must be in binary format and, in addition, can be Zlib compressed.
2. Operators can only be applied to VTU files sharing the same topology. The user must ensure that this condition is met.
3. The data type of numerical fields of the VTU file, for which operators should be applied, have to be `Float64`.

## Usage

The VTUFileHandler implements a basic VTU reader and writer through the functions:
```julia
function VTUFile(file::String) ... end 
function Base.write(vtu::VTUFile, add_timestamp=true) ... end
```
By default, a timestamp is added if VTU files are written to disk not to overwrite existing files. Only data fields that are registered by the function 
```julia
function set_uncompress_keywords(uk::Vector{String}) ... end
```
before reading the VTU file are uncompressed and can be altered. For applying math operators onto a data field, the associated field has to be registered by the function 
```julia
function set_interpolation_keywords(ik::Vector{String}) ... end
```
The following math operators are implemented:
```julia 
+(::VTUFile, ::VTUFile),+(::VTUFile, ::Number),
-(::VTUFile, ::VTUFile),-(::VTUFile, ::Number),
*(::VTUFile, ::VTUFile),*(::VTUFile, ::Number),
/(::VTUFile, ::VTUFile),/(::VTUFile, ::Number),
^(::VTUFile, ::Number)
```
In-place variations of the operators above are implemented as well.


## Example

A three-dimensional cube with dimension (x,y,z) with 0<=x,y,z<=2 discretized by quadrilian elements with 27 points and 8 cells named [`vox8.vtu`](https://github.com/baxmittens/VTUFileHandler.jl/blob/main/test/vox8.vtu) with a linear ramp in x-direction (f(x=0,y,z)=0, f(x=2,y,z)=0.8) as a result field termed `xramp` will be used as an example. The following set of instructions transform the result field from a linear ramp to a quadratic function in x-direction (displayed as a piecewise linear field due to the discretization):
```julia
using VTUFileHandler
set_uncompress_keywords(["xRamp"]) # uncrompress data field xramp
set_interpolation_keywords(["xRamp"]) # apply math operators to xramp
vtu = VTUFile("vox8.vtu"); # read the vtu
vtu += vtu/4; # [0.0,...,0.8] -> [0.0,...,1.0]
vtu *= 4.0; # [0,...,1.0] -> [0.0,...,4.0]
vtu -= 2.0; # [0,...,4.0] -> [-2.0,...,2.0]
vtu ^= 2.0; # [-2.0,...,2.0] -> [4.0,...,0.0,...,4.0]
rename!(vtu,"vox8_1.vtu")
write(vtu)
```

## Contributions, report bugs and support

Contributions to or questions about this project are welcome. Feel free to create a issue or a pull request on [GitHub](https://github.com/baxmittens/VTUFileHandler.jl).

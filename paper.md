---
title: 'DistributedSparseGrids.jl: A Julia library implementing an Adaptive Sparse Grid collocation method'
tags:
  - Julia
  - stochastics
  - sparse grids
  - high-performance computing
authors:
  - name: Maximilian Bittens
    orcid: 0000-0001-9954-294X
#   equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Robert L. Gates
#   equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: Federal Institute for Geosciences and Natural Resources (BGR)
   index: 1
 - name: Zeiss
   index: 2
date: 11 October 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Abstract

Computational integration or interpolation of high-dimensional functions is subject to the curse of dimensionality on full tensor grids. Sparse-Grid approximation can mitigate the latter, especially if the evaluation of the underlying function is costly. In the following, a Julia implementation of an local-lagrangian adaptive hierarchical sparse grid collocation method is presented, which is suitable for memory-heavy objects generated on distributed workers.

# Statement of need

[DistributedSparseGrids.jl](https://github.com/baxmittens/DistributedSparseGrids.jl) is a Julia package for integration and interpolation of functions with generic return types. There are other approaches to sparse grid approximation written in the julia language, as [SparseGrids.jl](https://github.com/robertdj/SparseGrids.jl), [AdaptiveSparseGrids.jl](https://github.com/jacobadenbaum/AdaptiveSparseGrids.jl), [GalerkinSparseGrids.jl](https://github.com/ABAtanasov/GalerkinSparseGrids.jl) or [Tasmanian.jl](https://github.com/floswald/Tasmanian.jl). However, there is no Julia package available at the moment which is suitable if the solution of the underlying (discretized) physical problem is time and resource consuming, requiring it to be solve on either a server or cluster enviroment.

# Introduction

Sparse tensor product quadrature rules, mitigating the curse of dimensionality occurring in full tensor grid constructions, were provided by `@smolyak1963quadrature`. In the last two decades collocation methods were prominent in the solution of
stochastic partial differential equation as shown in `@babuvska2007stochastic` and `@nobile2008sparse`.
`@ma2009adaptive` were able to once again increase efficiency of the collocation approach
by introducing an error-adaptive formulation of the method, which will serve as a basis for the
collocation method described in this project. For more information about the theory of the method implemented, see `@gates2015multilevel`.

# Features


- Nested one-dimensional Clenshaw-Curtis rule
- Smolyak's sparse grid construction
- local hierarchical Lagrangian basis
- different pointsets (open, closed, halfopen)
- adaptive refinement
- distributed function evaluation with ```Distributed.remotecall_fetch```
- multi-threaded calculation of basis coefficients with ```Threads.@threads```
- usage of arbitrary return types 
- integration
- experimental: integration over $X_{\sim (i)}$ (the $X_{\sim (i)}$  notation indicates the set of all variables except $X_{i}$).


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
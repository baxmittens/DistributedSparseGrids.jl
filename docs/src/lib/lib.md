# Library

## Contents 

```@contents
Pages = ["lib.md"]
Depth = 4
```

## Functions

### Index

```@index
Pages = ["lib.md"]
```

### Typedefs

```@docs
DistributedSparseGrids.AbstractSparseGrid
DistributedSparseGrids.AbstractHierarchicalSparseGrid
DistributedSparseGrids.PointDict
```

### Structs

```@docs
DistributedSparseGrids.CollocationPoint
DistributedSparseGrids.HierarchicalCollocationPoint
DistributedSparseGrids.AdaptiveHierarchicalSparseGrid
```

### General functions


```@docs
init
interpolate
integrate
init_weights!
init_weights_inplace_ops!
distributed_init_weights!
distributed_init_weights_inplace_ops!
generate_next_level!
```

### Utils

```@docs
DistributedSparseGrids.AHSG
```



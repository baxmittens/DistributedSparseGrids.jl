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

### General functions

```@docs
VTUHeader
VTUFileHandler.VTUDataField
VTUFileHandler.VTUData
VTUFile
```

### VTU Keywords

```@docs
VTUFileHandler.VTUKeyWords
set_uncompress_keywords
add_uncompress_keywords
set_interpolation_keywords
add_interpolation_keywords
```

### Utils

```@docs
getindex(::VTUFile, ::String)
VTUFileHandler.rename!
deepcopy(::VTUFile)
similar(::VTUFile)
fill!(::VTUFile, ::Float64)
zero(::VTUFile)
one(::VTUFile)
empty!(::VTUFile)
empty(::VTUFile)
```

### IO-Functions

```@docs
write(::VTUFile, ::Bool)
```
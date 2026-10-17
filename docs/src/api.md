# API

```@index
```

```@autodocs
Modules = [BICsCom]
```

## Internals

These functions are not exported and are not part of the public API;
they may change without notice.

```@docs
BICsCom._fix_phase!
BICsCom._reconstruct_dirichlet
BICsCom._edge_fourier_coeffs
```

## Extension: CMakieExt

The plotting functions 
- [`plot_square!`](@ref)
- [`plot_samplingpoints!`](@ref)
- [`plot_min_svals!`](@ref)
- [`plot_field_mode!`](@ref)
are implemented by the `CMakieExt` extension and become available once `CairoMakie`
is loaded. Their documentation is listed above with the rest of the package API.
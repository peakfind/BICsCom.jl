```@meta
CurrentModule = BICsCom
```

# Semi-analytic method

## DtN matrix based on Spherical expansion

We follow [Huang2006](@cite)

## Imposing boundary conditions

## Phase convention of the computed field

The eigenvector obtained by [`compute_mode`](@ref) is determined only up to a
complex factor, so the field returned by [`evaluate_field`](@ref) inherits an
arbitrary normalization and phase — one should not rely on the phase chosen by
the eigensolver. The normalization is fixed by [`normalize_field!`](@ref), and
the phase must be fixed by a physical convention before comparing the field
with published plots. For a propagating BIC, the natural choice is the
symmetry gauge: the phase (unique up to a sign) such that
``u(-x, y) = \overline{u(x, y)}`` holds for an x-even mode, with the remaining
sign fixed by requiring the dominant entry to be positive real. For standing
BICs this is equivalent to choosing the phase that makes the field real.
Plots of ``|u|`` are unaffected by the phase.

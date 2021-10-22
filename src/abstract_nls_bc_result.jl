# ================================================================
# Result abstraction for bound constrained problems
# ================================================================
#
# For the moment only use an alias as we do not have specific methods
#
@doc raw"""
```julia
const AbstractNLSBCResult = AbstractNLSResult
```

NLS with bound contraints solver result abstraction. 

For the moment only an **alias** of [`AbstractNLSResult`](@ref), with
the **same interface**.
"""
const AbstractNLSBCResult = AbstractNLSResult


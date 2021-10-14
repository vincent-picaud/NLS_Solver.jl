"""
```julia
@enum(BoundConstraint_Enum,
      ACTIVE_LB = -1,
      INACTIVE_BC = 0,
      ACTIVE_UB = +1)
```

An enum to store bound constraint state.
"""
@enum(BoundConstraint_Enum,
      ACTIVE_LB = -1,
      INACTIVE_BC = 0,
      ACTIVE_UB = +1)

# Regularization Schedule

A structure to store a regularization factor sequence. 

There is a base abstract class [`AbstractRegularizationSchedule`](@ref) with its interface:
- [`burning_phase`](@ref)
- [`regularization_factor`](@ref)


Two concrete implementations are provided:
- [`NoRegularizationSchedule`](@ref)
- [`ExpRegularizationSchedule`](@ref)

## Documentation

```@autodocs
Modules = [NLS_Solver]
Pages = ["regularization_schedule.jl"]
Private = false
```

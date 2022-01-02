```@meta
CurrentModule = NLS_Solver
```

```@setup session
using NLS_Solver
using Random
using Plots
ENV["GKSwstype"]=100
gr()
```

# Nonlinear regressions

This **NLS_Solver.jl** package is a generic nonlinear least squares
solver. However a common use case of this kind of solver is for
performing nonlinear regressions. In this tutorial we show how this
package can be used in this context. You can reproduce the computation
thanks to `sandbox/nonlinear_regression.jl`.

## Synthetic data

We will fit an exponential baseline plus a gaussian peak of the form:
```math
h \exp{\left(-\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2\right)} + b \exp{\left(-\frac{x}{\tau}\right)}
```

We first define some bricks for our model

```@example session
function gaussian_peak(x::Real,
                       h::Real,
                       μ::Real,
                       σ::Real)
    @assert σ>0
    
    h*exp(-((x-μ)/σ)^2/2)
end

function exp_baseline(x::Real,
                      b::Real,
                      τ::Real)
    @assert τ>0

    b*exp(-x/τ)
end
          
function model(x::Real,θ::AbstractVector)
    @assert length(θ) == 5
    
    # one gaussian peaks + exp baseline
    #
    (h, μ, σ, b, τ) = (θ_i for θ_i in θ)
    
    gaussian_peak(x,h,μ,σ) + exp_baseline(x,b,τ) 
end

nothing # hide
```

Now the important function that computes the residue:
```math
r_i(θ) = Y_i - m(X_i,θ),\ i=1,\dots,n_{\text{sample}}
```
where the ``r_i`` components are used to define the objective function ``\frac{1}{2}\|r(θ)\|_2^2``.

```@example session
function residue(X::AbstractVector,Y::AbstractVector,θ::AbstractVector)
    map(zip(X,Y)) do (X_i,Y_i)
        Y_i - model(X_i, θ)
    end
end

nothing # hide
```

We are now ready to generate our synthetic data:

```@example session
X=[1:0.1:30;]
θ_true=[2.0,15,3.0,3.0,8.0]

Y_true  = map(X_i -> model(X_i,θ_true), X)
Y_noisy = map(Y_i -> Y_i + 0.5 * randn(), Y_true)

plot(X,Y_noisy,label="noisy signal")
plot!(X,Y_true,linewidth=3,linestyle=:dot,label="true signal")
```

## Solve the problem

Now we solve the problem as we have done for [Bound constrained
nonlinear least squares](@ref bc_nls_pb).

The problem dimensions are:
```@example session
n_θ = length(θ_true)
n_sample = length(X)
nothing # hide
```

We impose some bound constraints:

```@example session
#                         h,  μ,  σ,   b,  τ
θ_lowerbound = Float64[   1, 10,  1,   0,  1 ]
θ_upperbound = Float64[ Inf, 20, 10, Inf, 10 ]
bc = BoundConstraints(θ_lowerbound, θ_upperbound)
nothing # hide
```

The initial guess is the  vector ``\mathbf{1}``. 

```@example session
θ_init = ones(n_θ)
nothing # hide
```

!!! note
    You can call the `solve()` function with an unfeasible initial vector θ. 
	(where =unfeasible= means that the bound constraints are violated). 

Then we choose the solver: 
```@example session
conf = LevenbergMarquardt_BC_Conf()
nothing # hide
```

We define the nonlinear least squares problem from the `residue()` function:
```@example session
nls = create_NLS_problem_using_ForwardDiff(θ->residue(X,Y_noisy,θ),n_θ => n_sample)
nothing # hide
```
The residue Jacobian is computed using automatic differentiation.

We can now solve the problem:

```@example session
result = solve(nls,θ_init,bc,conf)
nothing # hide
```

## Fit result

We check that the solver converged:

```@example session
@assert converged(result)
```

```@example session
θ_solution = solution(result)
```

and plot the fitted model:

```@example session
Y_solution = map(X_i -> model(X_i,θ_solution), X)
plot!(X,Y_solution,linewidth=3,label="fitted model")
```

We can get further information from the `result` structure:

```@example session
iteration_count(result)
```

```@example session
objective_value(result)
```

!!! note 
    In the future I will add gradient, Hessian and multipliers at
    the solution.  However, please note that in the constrained case
    you cannot use ``(J^tJ)^{-1}`` directly to estimate your
    parameter uncertainties.

## More fit examples

This **NLS_Solver.jl** package is used by the
[NLS_Fit.jl](https://github.com/vincent-picaud/NLS_Fit.jl) package
which is dedicated to peak fitting in spectrometry. You will find
there other examples of model fittings.

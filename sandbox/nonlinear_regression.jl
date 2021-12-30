# This example is reused in the documentation
# -> NLS_Solver.jl/docs/src/nonlinear_regression.md
#
using NLS_Solver
using Random
using Plots

# Be reproducible
rng = MersenneTwister(1234)

# Define model ================
#
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
    
    exp_baseline(x,b,τ) + gaussian_peak(x,h,μ,σ)
end

function residue(X::AbstractVector,Y::AbstractVector,θ::AbstractVector)
    map(zip(X,Y)) do (X_i,Y_i)
        Y_i - model(X_i, θ)
    end
end

# Generate synthetic data ================
#
X=[1:0.1:30;]
θ_true=[2.0,15,3.0,3.0,8.0]

Y_true = map(X) do X_i
    model(X_i,θ_true) 
end
Y_noisy = map(Y_true) do Y_i
    Y_i + 0.5 * randn()
end

plot(X,Y_noisy,label="noisy signal")
plot!(X,Y_true,linewidth=3,linestyle=:dot,label="true signal")

# Solve the problem ================
#
n_θ = length(θ_true)
n_sample = length(X)

# Bound constraints ----------------
#
#                         h,  μ,  σ,   b,  τ
θ_lowerbound = Float64[   1, 10,  1,   0,  1 ]
θ_upperbound = Float64[ Inf, 20, 10, Inf, 10 ]
bc = BoundConstraints(θ_lowerbound, θ_upperbound)

# Initial guess
#
θ_init = project!(ones(n_θ),bc)

# Solver ----------------
#
conf = LevenbergMarquardt_BC_Conf()

# Wrap residue ----------------
#
nls = create_NLS_problem_using_ForwardDiff(θ->residue(X,Y_noisy,θ),n_θ => n_sample)

# Solve it! ----------------
#
result = solve(nls,θ_init,bc,conf)

# Plot fitted model ================
#

# Initial model ----------------
#
Y_init = map(X_i -> model(X_i,θ_init), X)
plot!(X,Y_init,linewidth=3,label="initial model")


# Fitted model ----------------
#
θ_solution = solution(result)

Y_solution = map(X_i -> model(X_i,θ_solution), X)
plot!(X,Y_solution,linewidth=3,label="fitted model")



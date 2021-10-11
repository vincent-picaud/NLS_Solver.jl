export KunischRendlConf

# Configuration
#

"""

Kunisch-Rendl solver configuration

# Regularization 

The Kunisch method is quite sensitive to ill-conditioned Q
matrix

To tackle such kind of problems we introduce a decreasing regularization
enforcing diagonal dominance

`c` is the diagonal amplification factor >= 1 applied to the
diagonal. `c` is the value at the first iteration. `c=1` at iteration
`n≥k` and the amplification has no more effect. The decreasing
schedule is empirical and defined as follows:

    cₖ = exp(-log(c₀) TODO to fix

Note: convergence tests are performed only after `k` iterations, hence
      the algorithm perform at least `k` iterations

""" 
struct KunischRendlConf <: AbstractQuadSolverConf
    iter_max::Int
    
    # Regularization
    c::Float64
    k::Int

    function KunischRendlConf(;
                              iter_max=50,
                              c=1,
                              k=1)
        @assert c≥1
        @assert iter_max>k≥1
        new(Int(iter_max),Float64(c),Int(k))
    end
end

# Output 
#
struct KunischRendlResult <: AbstractQuadSolverResult
    converged::Bool
    iter::Int
    X_solution::AbstractVector
    τ_solution::AbstractVector
end 

# ****************************************************************
# Algorithm
# ****************************************************************
#

# Modify
#
# Q.x=q
# 
# x[i]=a if Z[i]=active_lb
# x[i]=b if Z[i]=active_ub
#
# into an uncontrained system
#
# \tilden  Q'.x=q'

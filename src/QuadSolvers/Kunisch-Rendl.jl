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
    c::Float64
    k::Int

    function KunischRendlConf(;c=1,k=1)
        @assert c≥1
        @assert k≥1
        new(Float64(c),Int(k))
    end
end

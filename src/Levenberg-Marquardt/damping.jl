# ================================================================
# Damping factor
# ================================================================
#
struct LM_Damping
    _μ::Float64
    _ν::Float64
    
    function LM_Damping(μ,ν)
        @assert μ≥0 "μ≥0  μ=$μ"
        @assert ν>0
        new(μ,ν)
    end

    function LM_Damping(μ)
        LM_Damping(μ,2)
    end
end

get_μ(d::LM_Damping) = d._μ

increase_μ(d::LM_Damping) = LM_Damping(d._μ * d._ν, 2* d._ν)

function update_μ(d::LM_Damping,ρ::Real)
    if ρ>0
        return LM_Damping(d._μ *max(1/3,1-(2ρ-1)^3))
    end
    
    increase_μ(d)
end

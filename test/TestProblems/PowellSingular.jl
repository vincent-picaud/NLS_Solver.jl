# PowellSingular function
#
# source: Nielsen, Hans Bruun and others, Damping parameter in Marquardt's method, 1999, IMM
#
function PowellSingular()
NLS_ForwardDiff(4,4) do θ
        T = eltype(θ)

        T[θ[1]+10*θ[2],
          sqrt(5)*(θ[3]-θ[4]),
          (θ[2]-2*θ[3])^2,
          sqrt(10)*(θ[1]-θ[4])^2]

    end
end


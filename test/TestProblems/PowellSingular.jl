# PowellSingular function
#
# source: Nielsen, Hans Bruun and others, Damping parameter in Marquardt's method, 1999, IMM
#
function PowellSingular()
    create_nls_using_forwarddiff(4,4) do θ
        SVector{4,eltype(θ)}(
            θ[1]+10*θ[2],
            sqrt(5)*(θ[3]-θ[4]),
            (θ[2]-2*θ[3])^2,
            sqrt(10)*(θ[1]-θ[4])^2)
    end
end


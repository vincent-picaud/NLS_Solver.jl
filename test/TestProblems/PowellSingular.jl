# PowellSingular function
#
# source: Nielsen, Hans Bruun and others, Damping parameter in Marquardt's method, 1999, IMM
#
function PowellSingular()
    create_nls_using_forwarddiff(Float64,4,4) do r,θ
        r[1] = θ[1]+10*θ[2]
        r[2] = sqrt(5)*(θ[3]-θ[4])
        r[3] = (θ[2]-2*θ[3])^2
        r[4] = sqrt(10)*(θ[1]-θ[4])^2
    end
end


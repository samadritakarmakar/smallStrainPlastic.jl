function createVoigtElasticTensor(E::Float64, ν::Float64)::Array{Float64, 2}
    c::Float64 = E/((1+ν)*(1-2*ν))
    C::Array{Float64, 2} = zeros(6,6)
    C = [1-ν ν ν 0 0 0;
         ν 1-ν ν 0 0 0;
         ν ν 1-ν 0 0 0;
         0 0 0 (1-2*ν)/2 0 0;
         0 0 0 0 (1-2*ν)/2 0;
         0 0 0 0 0 (1-2*ν)/2]
    C = c*C
    return C
end

function ProjectionTensor4(i::Int64, j::Int64, k::Int64, l::Int64)
    return (i == k ? 1.0 : 0.0)*(j == l ? 1.0 : 0.0)-1.0/3.0*(i == j ? 1.0 : 0.0)*(k == l ? 1.0 : 0.0)
end

function getProjectionTensor4()
    return SymmetricTensor{4,3,Float64}(ProjectionTensor4)
end

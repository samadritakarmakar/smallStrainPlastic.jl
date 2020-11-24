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

function get_σₘ_𝐬(σ_voigt)
    σ::SymmetricTensor{2,3, Float64, 6} = deepcopy(Tensors.fromvoigt(SymmetricTensor{2,3}, σ_voigt))
    σₘ::Float64 = tr(σ)
    #Deviatoric Stress
    σ -= 1.0/3.0*σₘ*one(SymmetricTensor{2,3})
    𝐬::Float64 = sqrt(3.0/2.0)*norm(σ)
    return σₘ, 𝐬
end

function get_ϵₘ_𝒆(ϵ)
    𝜺::SymmetricTensor{2,3, Float64, 6} = deepcopy(Tensors.fromvoigt(SymmetricTensor{2,3}, ϵ))
    ϵₘ::Float64 = tr(𝜺)
    #Deviatoric Stress
    𝜺 -= 1.0/3.0*ϵₘ*one(SymmetricTensor{2,3})
    𝒆::Float64 = sqrt(2.0/3.0)*norm(𝜺)
    return ϵₘ, 𝒆
end

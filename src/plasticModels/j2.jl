
function 𝒇_j2(σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)
    σ_y::Float64 = params.f
    σ::SymmetricTensor{2,3, Float64, 6} = Tensors.fromvoigt(SymmetricTensor{2,3}, σ_voigt)
    #Deviatoric Stress
    σ -= 1.0/3.0*tr(σ)*one(SymmetricTensor{2,3})
    f::Float64 = sqrt(3/2).*norm(σ)-(σ_y-q[1])
    return f #<= 0.0 ? 0.0 : f
end

function ∂𝒇_∂𝛔_j2!(∂f_∂σ::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)
    σ_y::Float64 = params.∂f_∂σ
    σ::SymmetricTensor{2,3, Float64, 6} = Tensors.fromvoigt(SymmetricTensor{2,3}, σ_voigt)
    #Deviatoric Stress
    σ -= 1.0/3.0*tr(σ)*one(SymmetricTensor{2,3})
    ∂f∂σ::SymmetricTensor{2,3, Float64, 6} = sqrt(1.5)*(1.0/norm(σ)*Tensors.dcontract(getProjectionTensor4(),σ))
    ∂f_∂σ .= Tensors.tovoigt(∂f∂σ)
    #println("∂f_∂σ = ",∂f_∂σ)
    return ∂f_∂σ
end

function ∂𝒇_∂𝐪_j2!(∂f_∂q::Array{Float64, 1}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1},  params::ModelParams)
    f::Float64 = 𝒇_j2(σ_voigt, q, params)
    ∂f_∂q[1,1] = f <= 0.0 ? 0.0 : 1.0
    return ∂f_∂q
end

𝚯_j2! = ∂𝒇_∂𝛔_j2!

function ∂𝚯_∂𝛔_j2!(∂Θ_∂σ::Array{Float64, 2}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1}, params::ModelParams)
    ∂f_∂σ::Array{Float64, 1} = zeros(6)
    func(σ_voigt) = ∂𝒇_∂𝛔_j2!(∂f_∂σ, σ_voigt, q, params)
    denseJacobian!(∂Θ_∂σ, func, σ_voigt)
    return ∂Θ_∂σ
end

function ∂𝚯_∂𝐪_j2!(∂Θ_∂q::Array{Float64, 2}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1}, params::ModelParams)
    ∂Θ_∂q[1,1] = 0.0
    return ∂Θ_∂q
end

function 𝐡_j2!(h::Array{Float64, 1}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1}, params::ModelParams)

    h[1] = 0.0
    return h
end

function ∂𝐡_∂𝛔_j2!(∂h_∂σ::Array{Float64, 2}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1}, params::ModelParams)

    ∂h_∂σ[1,1] = 0.0
    return ∂h_∂σ
end

function ∂𝐡_∂𝐪_j2!(∂h_∂q::Array{Float64, 2}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1},  params::ModelParams)

    ∂h_∂q[1,1] = 0.0
    return ∂h_∂q
end

function 𝓗_j2!(H::Array{Float64, 1}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1}, α::Array{Float64, 1}, params::ModelParams)
    ∂f_∂q::Array{Float64, 1} = zeros(1)
    ∂𝒇_∂𝐪_j2!(∂f_∂q, σ_voigt, q, params)
    H[1] = params.H*∂f_∂q[1]
    return H
end

function ℂ_j2!(C::Array{Float64,2}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1},  params::ModelParams)

    ##Do nothing
    #C[:,:] = kwargs_C[1]
    return C
end

function 𝔻_j2!(D::Array{Float64,2}, σ_voigt::Array{Float64, 1},
    q::Array{Float64, 1},  params::ModelParams)

    D[1,1] = 1.0
    return D
end

j2Model = PlasticModel(𝒇_j2, ∂𝒇_∂𝛔_j2!, ∂𝒇_∂𝐪_j2!, 𝚯_j2!,
∂𝚯_∂𝛔_j2!, ∂𝚯_∂𝐪_j2!, 𝐡_j2!, ∂𝐡_∂𝛔_j2!, ∂𝐡_∂𝐪_j2!,
𝓗_j2!, ℂ_j2!, 𝔻_j2!, 6, 1)

function initParams_j2(σ_y::Float64,  params_H::Float64)
    params_f::Float64 = σ_y
    params_∂f_∂σ::Float64 = σ_y
    return ModelParams{Float64, Float64, Int64, Int64,
        Int64, Int64, Int64, Int64, Int64,
        Float64, Int64, Int64}(params_f, params_∂f_∂σ, 0,0,0,0,
        0, 0,0, params_H,0,0)
end

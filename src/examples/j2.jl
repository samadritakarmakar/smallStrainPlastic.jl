using SmallStrainPlastic, Tensors

function 𝒇_j2(σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}; kwargs...) where {T}
    σ_y::Float64 = kwargs[1]
    σ::SymmetricTensor{2,3, Float64, 6} = Tensors.fromvoigt(SymmetricTensor{2,3}, σ_voigt)
    #Deviatoric Stress
    σ -= 1.0/3.0*tr(σ)*one(SymmetricTensor{2,3})
    return sqrt(3/2).*norm(σ)-(σ_y-q[1])
end

function ∂𝒇_∂𝛔_j2!(∂f_∂σ, σ_voigt, q::Array{Float64, 1}; kwargs...)
    f(σ_voigt) = 𝒇_j2(σ_voigt, q; σ_y = kwargs[1])
    SmallStrainPlastic.denseJacobian!(∂f_∂σ, f, σ_voigt)
    return nothing
end

function ∂𝒇_∂𝐪_j2!(∂f_∂q::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}; kwargs...)
    ∂f_∂q[1,1] = 1.0
    return nothing
end

∂Θ_∂𝛔_j2! = ∂𝒇_∂𝛔_j2!
∂Θ_∂𝐪_j2! = ∂𝒇_∂𝐪_j2!

function ∂𝐡_∂𝛔_j2!(∂h_∂σ::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}; kwargs...)
    ∂h_∂σ[1,1] = 0.0
    return nothing
end

function ∂𝐡_∂𝐪_j2!(∂h_∂σ::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}; kwargs...)
    ∂h_∂σ[1,1] = 0.0
    return nothing
end

function 𝔻_j2!(D::Array{Float64,2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}; kwargs...)
    D[1,1] = 0.0
    return nothing
end

function ℂ_j2!(C::Array{Float64,2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}; kwargs...)
    C[:,:] = kwargs[1]
    return nothing
end

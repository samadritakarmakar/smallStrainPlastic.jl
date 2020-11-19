"""A plasticity new plasticity model can be defined by defining the following functions:

Yield Function:-

    𝒇(σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)

The partial of Yield Function with respect to stress, ∂𝒇/∂𝛔:-

    ∂𝒇_∂𝛔!(∂f_∂σ::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)

The partial of Yield Function with respect to hardening variable, ∂𝒇/∂𝐪:

    ∂𝒇_∂𝐪!(∂f_∂q::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1},  params::ModelParams)

The function in the flow rule for the plastic strain, 𝛆̇ᵖ = λ̇  𝚯(𝛔, 𝐪) :-

    𝚯!(Θ::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)

The partial of plastic strain flow rule function with respect to stress, ∂𝚯/∂𝛔:-

    ∂𝚯_∂𝛔!(∂Θ_∂σ::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)

The partial of plastic strain flow rule function with respect to hardening, ∂𝚯/∂𝐪:-

    ∂𝚯_∂𝐪!(∂Θ_∂q::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)

The function in the flow rule for the internal variable 𝛂̇ = λ̇  𝐡(𝛔, 𝐪):-

    𝐡!(h::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)

The partial of plastic strain flow rule function with respect to stress, ∂𝐡/∂𝛔:-

    ∂𝐡_∂𝛔!(∂h_∂σ::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)

The partial of plastic strain flow rule function with respect to hardening, ∂𝐡/∂𝐪:-

    ∂𝐡_∂𝐪!(∂h_∂q::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1},  params::ModelParams)

If the evolution of the hardening variable 𝐪̇ is defined as 𝐪̇ = -𝓗(𝛂), then the function it is
dependent on can be written as:-

    𝓗!(H::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, α::Array{Float64, 1}, params::ModelParams)

For ease of use, defining a function that saves the stiffness tensor is also made available :-

    ℂ!(C::Array{Float64,2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1},  params::ModelParams)

If the hardening variable 𝐪̇ is defined as 𝐪̇ = -𝓗(𝛂), then an equivalent to stiffness tensor defined as
ℂ = ∂𝛔/∂𝛆ᵉ, we can defined as 𝔻 = -∂𝐪/∂𝛂 = ∂𝓗(𝛂)/∂𝐪 :-

    𝔻!(D::Array{Float64,2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1},  params::ModelParams)
"""
struct PlasticModel
    𝒇::Function
    ∂𝒇_∂𝛔!::Function
    ∂𝒇_∂𝐪!::Function
    𝚯!::Function
    ∂𝚯_∂𝛔!::Function
    ∂𝚯_∂𝐪!::Function
    𝐡!::Function
    ∂𝐡_∂𝛔!::Function
    ∂𝐡_∂𝐪!::Function
    𝓗!::Function
    ℂ!::Function
    𝔻!::Function
    ϵSize::Int64
    αSize::Int64
end
"""This is a list of Plastic Variables most commonly used in the code.
The idea is to reduce the number of parameters passed to a function

    PlasticVars(C, σ_voigt, ϵ, ϵᵖ, H, q, α, Cᵀ)
    """
mutable struct PlasticVars
    C::Array{Float64, 2}
    D::Array{Float64, 2}
    σ_voigt::Array{Float64, 1}
    ϵ::Array{Float64, 1}
    ϵᵖ::Array{Float64, 1}
    H::Array{Float64, 1}
    q::Array{Float64, 1}
    α::Array{Float64, 1}
    Cᵀ::Array{Float64, 2}
end

"This function initializes Plastic variables. Since the number of
variables in Plasticity are quite a few, this function helps in
giving them initial values.

    plasticVars = initPlasticVars()
"
function initPlasticVars(model::PlasticModel)
    ϵSize::Int64 = model.ϵSize
    αSize::Int64 = model.αSize
    C::Array{Float64, 2} = zeros(ϵSize, ϵSize)
    D::Array{Float64, 2} = zeros(αSize, αSize)
    σ_voigt::Array{Float64, 1} = zeros(ϵSize)
    ϵ::Array{Float64, 1} = zeros(ϵSize)
    ϵᵖ::Array{Float64, 1} = zeros(ϵSize)
    H::Array{Float64, 1} = zeros(αSize)
    q::Array{Float64, 1} = zeros(αSize)
    α::Array{Float64, 1} = zeros(αSize)
    Cᵀ::Array{Float64, 2} = zeros(ϵSize, ϵSize)
    return PlasticVars(C, D, σ_voigt, ϵ, ϵᵖ, H, q, α, Cᵀ)
end

"""Abstract type Parameters to make it easy to pass variables to functions"""
abstract type  Parameters end

"""This is one common structure providing for the possiblity of the various
parameters that may be there for different models. It is recommended that
    each model definition has it's own definition of Initilization of these
    parameters so that a huge definition of modelParams can be avoided."""
struct ModelParams{params_f, params_∂f_∂σ, params_∂f_∂q, params_Θ,
    params_∂Θ_∂σ, params_∂Θ_∂q, params_h, params_∂h_∂σ, params_∂h_∂q,
    params_H, params_C, params_D} <: Parameters
    f::params_f
    ∂f_∂σ::params_∂f_∂σ
    ∂f_∂q::params_∂f_∂q
    Θ::params_Θ
    ∂Θ_∂σ::params_∂Θ_∂σ
    ∂Θ_∂q::params_∂Θ_∂q
    h::params_h
    ∂h_∂σ::params_∂h_∂σ
    ∂h_∂q::params_∂h_∂q
    H::params_H
    C::params_C
    D::params_D
end

mutable struct tolerances
    f::Float64
    R::Float64
end

"""This function is a that checks if the stess state is in the plastic state or not.
If in the plastic state it initiates the return mapping algorithm so that the stress
state is brought back on to surface of the yield surface.

    checkPlasticState!(plasticVars, model, elementNo integrationPt, parameters)
"""
function checkPlasticState!(plasticVars::PlasticVars, model::PlasticModel,
    params::ModelParams, elementNo::Int64, integrationPt::Int64)

    getState!(plasticVars.ϵᵖ, plasticVars.α, stateDict, elementNo, integrationPt)
    model.ℂ!(plasticVars.C, plasticVars.σ_voigt, plasticVars.q, params)
    plasticVars.σ_voigt .= plasticVars.C*(plasticVars.ϵ .- plasticVars.ϵᵖ)
    model.𝓗!(plasticVars.H, plasticVars.σ_voigt, plasticVars.q, plasticVars.α, params)
    #for i ∈ 1:length(plasticVars.q)
    #    plasticVars.q[i] = -plasticVars.H[i]
    #end
    plasticVars.q .= -plasticVars.H
    if model.𝒇(plasticVars.σ_voigt, plasticVars.q, params) > 0
        println("In plastic regime")
        returnMapping!(plasticVars, model, params)
        updateStateDict!(plasticVars.ϵᵖ, plasticVars.α, stateDict, elementNo, integrationPt)
        return true
    else
        println("In elastic regime")
        plasticVars.Cᵀ = plasticVars.C
        return false
    end
end

function findResidual(R::Array{Float64, 1}, plasticVars::PlasticVars, model::PlasticModel,
    params::ModelParams)

end

function initReturnMappingVars(model::PlasticModel)
    ∂f_∂σ::Array{Float64, 1} = zeros(model.ϵSize)
    ∂f_∂q::Array{Float64, 1} = zeros(model.αSize)
    ∂Θ_∂σ::Array{Float64, 2} = zeros(model.ϵSize, model.ϵSize)
    ∂Θ_∂q::Array{Float64, 2} = zeros(model.ϵSize, model.αSize)
    ∂h_∂σ::Array{Float64, 2} = zeros(model.αSize, model.ϵSize)
    ∂h_∂q::Array{Float64, 2} = zeros(model.αSize, model.αSize)
    A::Array{Float64, 2} = zeros(model.ϵSize+model.αSize, model.ϵSize+model.αSize)
    Θ::Array{Float64, 1} = zeros(model.ϵSize)
    h::Array{Float64, 1} = zeros(model.αSize)
    R::Array{Float64, 1} = ones(model.ϵSize+model.αSize)
    f::Float64 = 1.0
    Δλ::Float64 = 0.0
    dΔλ::Float64 = 0.0
    return ∂f_∂σ, ∂f_∂q, ∂Θ_∂σ, ∂Θ_∂q, ∂h_∂σ, ∂h_∂q, Θ, h, R, A, f, Δλ, dΔλ
end

function updateReturnMappingVars!(∂f_∂σ::Array{Float64, 1},
    ∂f_∂q::Array{Float64, 1},
    ∂Θ_∂σ::Array{Float64, 2}, ∂Θ_∂q::Array{Float64, 2},
    ∂h_∂σ::Array{Float64, 2}, ∂h_∂q::Array{Float64, 2},
    Θ::Array{Float64, 1}, h::Array{Float64, 1},
    plasticVars::PlasticVars, model::PlasticModel,
    params::ModelParams)

    model.∂𝒇_∂𝛔!(∂f_∂σ, plasticVars.σ_voigt, plasticVars.q, params)
    model.∂𝒇_∂𝐪!(∂f_∂q, plasticVars.σ_voigt, plasticVars.q, params)
    model.𝚯!(Θ, plasticVars.σ_voigt, plasticVars.q, params)
    model.∂𝚯_∂𝛔!(∂Θ_∂σ, plasticVars.σ_voigt, plasticVars.q, params)
    model.∂𝚯_∂𝐪!(∂Θ_∂q, plasticVars.σ_voigt, plasticVars.q, params)
    model.𝐡!(h, plasticVars.σ_voigt, plasticVars.q, params)
    model.𝓗!(plasticVars.H, plasticVars.σ_voigt, plasticVars.q, plasticVars.α, params)
    model.∂𝐡_∂𝛔!(∂h_∂σ, plasticVars.σ_voigt, plasticVars.q, params)
    model.∂𝐡_∂𝐪!(∂h_∂q, plasticVars.σ_voigt, plasticVars.q, params)
    model.ℂ!(plasticVars.C, plasticVars.σ_voigt, plasticVars.q, params)
    model.𝔻!(plasticVars.D, plasticVars.σ_voigt, plasticVars.q, params)
    return nothing
end

function returnMapping!(plasticVars::PlasticVars, model::PlasticModel,
    params::ModelParams)

    ∂f_∂σ::Array{Float64, 1}, ∂f_∂q::Array{Float64, 1},
    ∂Θ_∂σ::Array{Float64, 2}, ∂Θ_∂q::Array{Float64, 2},
    ∂h_∂σ::Array{Float64, 2}, ∂h_∂q::Array{Float64, 2},
    Θ::Array{Float64, 1}, h::Array{Float64, 1},
    R::Array{Float64, 1}, A::Array{Float64, 2},
    f::Float64, Δλ::Float64, dΔλ::Float64 = initReturnMappingVars(model)
    ϵᵖα_n1::Array{Float64, 1} = [plasticVars.ϵᵖ; plasticVars.α]
    while (norm(f)>1e-12)#|| norm(R)/norm(plasticVars.ϵ)>1e-7)

        #Update Return mapping internal arrays
        updateReturnMappingVars!(∂f_∂σ, ∂f_∂q, ∂Θ_∂σ, ∂Θ_∂q, ∂h_∂σ, ∂h_∂q, Θ, h, plasticVars, model, params)
        #Update Residual
        plasticVars.σ_voigt .= plasticVars.C*(plasticVars.ϵ - plasticVars.ϵᵖ)
        plasticVars.q .= -plasticVars.H
        f = model.𝒇(plasticVars.σ_voigt, plasticVars.q, params)
        Θh::Array{Float64, 1} = [Θ; h]
        R .= -ϵᵖα_n1 + [plasticVars.ϵᵖ; plasticVars.α] + Δλ*Θh
        #update matrix [A]
        A[1:model.ϵSize,1:model.ϵSize] .= inv(plasticVars.C) + Δλ*∂Θ_∂σ
        A[model.ϵSize+1:model.ϵSize+model.αSize, 1:model.ϵSize] .= Δλ*∂h_∂σ
        A[1:model.ϵSize, model.ϵSize+1:model.ϵSize+model.αSize] .= Δλ*∂Θ_∂q
        A[model.ϵSize+1:model.ϵSize+model.αSize, model.ϵSize+1:model.ϵSize+model.αSize] .=
        inv(plasticVars.D)+ Δλ*∂h_∂q

        fA = [∂f_∂σ..., ∂f_∂q...]'*A
        dΔλ = (f .- fA*R)/(fA*Θh)
        Δλ += dΔλ
        plasticVars.ϵᵖ .= ϵᵖα_n1[1:model.ϵSize]
        plasticVars.α .= ϵᵖα_n1[model.ϵSize+1:model.ϵSize+model.αSize]
        C_D_inv::Array{Float64, 2} = [inv(plasticVars.C) zeros(model.ϵSize, model.αSize); zeros(model.αSize, model.ϵSize) inv(plasticVars.D)]
        ϵᵖα_n1 +=C_D_inv*A*(R .+ dΔλ*Θh)
        println("f = ", f, " norm(R)/norm(ϵ) = ", norm(R)/norm(plasticVars.ϵ), " dΔλ = ", dΔλ)
    end
end

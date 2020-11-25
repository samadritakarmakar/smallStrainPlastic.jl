"""This function is a that checks if the stess state is in the plastic state or not.
If in the plastic state it initiates the return mapping algorithm so that the stress
state is brought back on to surface of the yield surface.

    checkPlasticState!(plasticVars, model, elementNo integrationPt, parameters)
"""
function checkPlasticState!(plasticVars::PlasticVars, model::PlasticModel,
    params::ModelParams, elementNo::Int64, integrationPt::Int64)

    getState!(plasticVars.ϵᵖ, plasticVars.α, stateDict, elementNo, integrationPt)
    model.ℂ!(plasticVars.C, plasticVars.σ_voigt, plasticVars.q, plasticVars, params)
    plasticVars.σ_voigt = plasticVars.C*(plasticVars.ϵ .- plasticVars.ϵᵖ)
    model.𝓗!(plasticVars.H, plasticVars.σ_voigt, plasticVars.q, plasticVars.α, plasticVars, params)
    #for i ∈ 1:length(plasticVars.q)
    #    plasticVars.q[i] = -plasticVars.H[i]
    #end
    plasticVars.q = -plasticVars.H
    if model.𝒇(plasticVars.σ_voigt, plasticVars.q, plasticVars, params) > 0
        #println("In plastic regime")
        returnMapping!(plasticVars, model, params)
        updateStateDict!(plasticVars.ϵᵖ, plasticVars.α, stateDictBuffer, elementNo, integrationPt)
        return true
    else
        #println("In elastic regime")
        plasticVars.Cᵀ = plasticVars.C
        return false
    end
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

    model.∂𝒇_∂𝛔!(∂f_∂σ, plasticVars.σ_voigt, plasticVars.q,  plasticVars, params)
    model.∂𝒇_∂𝐪!(∂f_∂q, plasticVars.σ_voigt, plasticVars.q,  plasticVars, params)
    model.𝚯!(Θ, plasticVars.σ_voigt, plasticVars.q,  plasticVars, params)
    model.∂𝚯_∂𝛔!(∂Θ_∂σ, plasticVars.σ_voigt, plasticVars.q,  plasticVars, params)
    model.∂𝚯_∂𝐪!(∂Θ_∂q, plasticVars.σ_voigt, plasticVars.q, plasticVars, params)
    model.𝐡!(h, plasticVars.σ_voigt, plasticVars.q, plasticVars, params)
    model.∂𝐡_∂𝛔!(∂h_∂σ, plasticVars.σ_voigt, plasticVars.q, plasticVars, params)
    model.∂𝐡_∂𝐪!(∂h_∂q, plasticVars.σ_voigt, plasticVars.q, plasticVars, params)
    model.ℂ!(plasticVars.C, plasticVars.σ_voigt, plasticVars.q, plasticVars, params)
    model.𝔻!(plasticVars.D, plasticVars.σ_voigt, plasticVars.q, plasticVars, params)
    return nothing
end

"""This function is responsible for executing the return mapping algorithm. It
does so by calculating the evolution of the plastic strain using the Closest Point Projection
method. The following formulations are used

Here's an equation:

``d(\\Delta\\lambda) = \\frac{f^k - \\begin{bmatrix}\\partial f^k/\\partial \\sigma & \\partial f^k/ \\partial q\\end{bmatrix}\\begin{bmatrix}A\\end{bmatrix}\\begin{bmatrix} R \\end{bmatrix}}
{\\begin{bmatrix}\\partial f^k/\\partial \\sigma & \\partial f^k/ \\partial q\\end{bmatrix}\\begin{bmatrix}A\\end{bmatrix}\\begin{bmatrix}\\Theta \\\\ h \\end{bmatrix}}``

where:

``\\begin{bmatrix}R\\end{bmatrix} = -\\begin{bmatrix} \\epsilon^p_{n+1} \\\\ \\alpha_{n+1} \\end{bmatrix}
+\\begin{bmatrix} \\epsilon^p_{n} \\\\ \\alpha_{n} \\end{bmatrix}
+\\Delta\\lambda\\begin{bmatrix} \\Theta(\\sigma_{n+1}, q_{n+1}) \\\\ h(\\sigma_{n+1}, q_{n+1}) \\end{bmatrix}``

``\\begin{bmatrix}A\\end{bmatrix}^{-1} =
\\begin{bmatrix} \\bm{C}^{-1} + \\Delta\\lambda \\frac{\\partial \\Theta}{\\partial\\sigma_{n+1}} &
\\Delta\\lambda \\frac{\\partial \\Theta}{\\partial q_{n+1}}
\\\\ \\Delta\\lambda \\frac{\\partial h}{\\partial\\sigma_{n+1}} &
\\bm{D}^{-1} + \\Delta\\lambda \\frac{\\partial h}{\\partial q_{n+1}}
\\end{bmatrix}``

The Strain 𝛆ᵖ and the internal variable 𝛂 are updated as,

``\\begin{bmatrix}\\Delta \\epsilon^p \\\\ \\Delta \\alpha \\end{bmatrix} =
\\begin{bmatrix}\\bm{C}^{-1} & 0 \\\\ 0 & \\bm{D}^{-1} \\end{bmatrix}
\\begin{bmatrix}A\\end{bmatrix}
\\begin{bmatrix}\\Theta \\\\ h \\end{bmatrix}
d(\\Delta\\lambda)``
"""
function returnMapping!(plasticVars::PlasticVars, model::PlasticModel,
    params::ModelParams)

    ∂f_∂σ::Array{Float64, 1}, ∂f_∂q::Array{Float64, 1},
    ∂Θ_∂σ::Array{Float64, 2}, ∂Θ_∂q::Array{Float64, 2},
    ∂h_∂σ::Array{Float64, 2}, ∂h_∂q::Array{Float64, 2},
    Θ::Array{Float64, 1}, h::Array{Float64, 1},
    R::Array{Float64, 1}, A::Array{Float64, 2},
    f::Float64, Δλ::Float64, dΔλ::Float64 = initReturnMappingVars(model)

    iter::Int64 = 0
    ϵᵖα_n1::Array{Float64, 1} = [plasticVars.ϵᵖ; plasticVars.α]
    while ((norm(f)> tolerance.f|| norm(R)> tolerance.R) && iter<=tolerance.maxIter)
        ϵₘ::Float64, plasticVars.𝒆ᵖ = get_ϵₘ_𝒆(ϵᵖα_n1[1:model.ϵSize])
        plasticVars.σ_voigt = plasticVars.C*(plasticVars.ϵ - ϵᵖα_n1[1:model.ϵSize])
        model.𝓗!(plasticVars.H, plasticVars.σ_voigt, plasticVars.q,
        ϵᵖα_n1[model.ϵSize+1:model.ϵSize+model.αSize], plasticVars, params)
        plasticVars.q = -plasticVars.H
        f = model.𝒇(plasticVars.σ_voigt, plasticVars.q, plasticVars, params)
        #Update Return mapping internal arrays
        updateReturnMappingVars!(∂f_∂σ, ∂f_∂q, ∂Θ_∂σ, ∂Θ_∂q, ∂h_∂σ, ∂h_∂q, Θ, h, plasticVars, model, params)
        Θh::Array{Float64, 1} = [Θ; h]
        #Update Residual
        R = -ϵᵖα_n1 + [plasticVars.ϵᵖ; plasticVars.α] + Δλ*Θh
        #println("-ϵᵖα_n1 ", -ϵᵖα_n1," Δλ*Θh = ", Δλ*Θh)
        #update matrix [A]
        A[1:model.ϵSize,1:model.ϵSize] = inv(plasticVars.C) + Δλ*∂Θ_∂σ
        A[model.ϵSize+1:model.ϵSize+model.αSize, 1:model.ϵSize] = Δλ*∂h_∂σ
        A[1:model.ϵSize, model.ϵSize+1:model.ϵSize+model.αSize] = Δλ*∂Θ_∂q
        A[model.ϵSize+1:model.ϵSize+model.αSize, model.ϵSize+1:model.ϵSize+model.αSize] =
        inv(plasticVars.D)+ Δλ*∂h_∂q

        A = inv(A)
        
        fA = [∂f_∂σ..., ∂f_∂q...]'*A
        dΔλ = (f - fA*R)/(fA*Θh)
        Δλ += dΔλ
        C_D_inv::Array{Float64, 2} = inv([(plasticVars.C) zeros(model.ϵSize, model.αSize);
                                        zeros(model.αSize, model.ϵSize) (plasticVars.D)])
        ϵᵖα_n1 +=C_D_inv*A*(R + dΔλ*Θh)
        iter += 1
        println("f = ", f, " norm(R) = ", norm(R), " dΔλ = ", dΔλ)
    end
    if iter > tolerance.maxIter
        @warn "Return Mapping Exited without convergence"
    end
    plasticVars.ϵᵖ = ϵᵖα_n1[1:model.ϵSize]
    plasticVars.α = ϵᵖα_n1[model.ϵSize+1:model.ϵSize+model.αSize]
    plasticVars.σ_voigt = plasticVars.C*(plasticVars.ϵ - plasticVars.ϵᵖ)
    model.𝓗!(plasticVars.H, plasticVars.σ_voigt, plasticVars.q,
    plasticVars.α, plasticVars, params)
end

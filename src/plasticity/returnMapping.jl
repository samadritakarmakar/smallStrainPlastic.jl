"""This function is a that checks if the stess state is in the plastic state or not.
If in the plastic state it initiates the return mapping algorithm so that the stress
state is brought back on to surface of the yield surface.
    checkPlasticState!(plasticVars, model, elementNo integrationPt, parameters)
"""
function checkPlasticState!(plasticVars::PlasticVars, model::PlasticModel,
    params::Parameters, stateDict::Dict{T}, stateDictBuffer::Dict{T},
     elementNo::Int64, integrationPt::Int64;
     tolerance::Tolerance = Tolerance(1e-8, 1e-8, 1000), algoTangent = false) where T

    getState!(plasticVars.ϵᵖ, plasticVars.α, stateDict, elementNo, integrationPt)
    plasticVars.C .= model.ℂ(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
    if length(plasticVars.ϵ) == 6
        plasticVars.ϵ = get_P()*plasticVars.ϵ
    end
    plasticVars.σ_mandel .= plasticVars.C*(plasticVars.ϵ .- plasticVars.ϵᵖ)
    plasticVars.H .= model.𝓗(plasticVars.σ_mandel, plasticVars.q, plasticVars.α, plasticVars, params)
    plasticVars.q .= -plasticVars.H
    σ = returnMapping!(plasticVars, model, params, tolerance = tolerance, algoTangent = algoTangent)
    updateStateDict!(plasticVars.ϵᵖ, plasticVars.α, stateDictBuffer,
    elementNo, integrationPt)
    return σ
end


#=function initReturnMappingVars(model::PlasticModel)
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
end=#

#=function updateReturnMappingVars!(∂f_∂σ::Array{Float64, 1},
    ∂f_∂q::Array{Float64, 1},
    ∂Θ_∂σ::Array{Float64, 2}, ∂Θ_∂q::Array{Float64, 2},
    ∂h_∂σ::Array{Float64, 2}, ∂h_∂q::Array{Float64, 2},
    Θ::Array{Float64, 1}, h::Array{Float64, 1},
    plasticVars::PlasticVars, model::PlasticModel,
    params::ModelParams)

    model.∂𝒇_∂𝛔!(∂f_∂σ, plasticVars.σ_mandel, plasticVars.q,  plasticVars, params)
    model.∂𝒇_∂𝐪!(∂f_∂q, plasticVars.σ_mandel, plasticVars.q,  plasticVars, params)
    model.𝚯!(Θ, plasticVars.σ_mandel, plasticVars.q,  plasticVars, params)
    model.∂𝚯_∂𝛔!(∂Θ_∂σ, plasticVars.σ_mandel, plasticVars.q,  plasticVars, params)
    model.∂𝚯_∂𝐪!(∂Θ_∂q, plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
    model.𝐡!(h, plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
    model.∂𝐡_∂𝛔!(∂h_∂σ, plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
    model.∂𝐡_∂𝐪!(∂h_∂q, plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
    model.ℂ!(plasticVars.C, plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
    model.𝔻!(plasticVars.D, plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
    return nothing
end=#

function calculateMatrix_A(Δλ::Float64, plasticVars::PlasticVars, model::PlasticModel,
    params::Parameters)
    A::Array{Float64, 2} = zeros(model.ϵSize+model.αSize, model.ϵSize+model.αSize)
    #println("plasticVars.C = ", plasticVars.C)
    A[1:model.ϵSize,1:model.ϵSize] = inv(plasticVars.C) +
    Δλ*model.∂𝚯_∂𝛔(plasticVars.σ_mandel, plasticVars.q,  plasticVars, params)

    A[model.ϵSize+1:model.ϵSize+model.αSize, 1:model.ϵSize] = Δλ*
    model.∂𝐡_∂𝛔(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)

    A[1:model.ϵSize, model.ϵSize+1:model.ϵSize+model.αSize] = Δλ*
    model.∂𝚯_∂𝐪(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)

    A[model.ϵSize+1:model.ϵSize+model.αSize, model.ϵSize+1:model.ϵSize+model.αSize] =
    inv(plasticVars.D)+ Δλ*model.∂𝐡_∂𝐪(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)

    return inv(A)
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
    params::Parameters; tolerance::Tolerance = Tolerance(1e-8, 1e-8, 1000), algoTangent = false)
    Δλ::Float64 = 0.0
    dΔλ::Float64 = 0.0
    f = model.𝒇(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
    if(f>0)
        #updateReturnMappingVars!(∂f_∂σ, ∂f_∂q, ∂Θ_∂σ, ∂Θ_∂q, ∂h_∂σ, ∂h_∂q, Θ, h, plasticVars, model, params)
        #Θh = [Θ; h]
        iter::Int64 = 0
        ϵᵖα_n1::Array{Float64, 1} = [plasticVars.ϵᵖ; plasticVars.α]
        Θh::Array{Float64, 1} = [model.𝚯(plasticVars.σ_mandel, plasticVars.q,  plasticVars, params);
        model.𝐡(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)]

        fA = [zeros(model.ϵSize)' zeros(model.αSize)']
        #Update Residual
        R = -ϵᵖα_n1 + [plasticVars.ϵᵖ; plasticVars.α] + Δλ*Θh
        while ((f> tolerance.f|| norm(R)> tolerance.R) && iter<=tolerance.maxIter)

            plasticVars.C .= model.ℂ(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
            plasticVars.D .= model.𝔻(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
            A = calculateMatrix_A(Δλ, plasticVars, model, params)
            fA .= [model.∂𝒇_∂𝛔(plasticVars.σ_mandel, plasticVars.q,  plasticVars, params);
            model.∂𝒇_∂𝐪(plasticVars.σ_mandel, plasticVars.q,  plasticVars, params)]'*A

            dΔλ = (f .- fA*R)/(fA*Θh)
            Δλ += dΔλ
            C_D_inv::Array{Float64, 2} = ([inv(plasticVars.C) zeros(model.ϵSize, model.αSize);
                                            zeros(model.αSize, model.ϵSize) inv(plasticVars.D)])
            Δσ_Δα = -A*(R + dΔλ*Θh)
            ϵᵖα_n1 += -C_D_inv*Δσ_Δα
            plasticVars.σ_mandel += Δσ_Δα[1:model.ϵSize]
            plasticVars.q += Δσ_Δα[model.ϵSize+1:model.ϵSize+model.αSize]
            f = model.𝒇(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
            Θh = [model.𝚯(plasticVars.σ_mandel, plasticVars.q,  plasticVars, params);
            model.𝐡(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)]

            #Update Residual
            R .= -ϵᵖα_n1 + [plasticVars.ϵᵖ; plasticVars.α] + Δλ*Θh
            iter += 1
            #println("f = ", f, " norm(R) = ", norm(R), " dΔλ = ", dΔλ)
        end
        if iter > tolerance.maxIter
            @warn "Return Mapping Exited without convergence"
        end

        plasticVars.ϵᵖ = ϵᵖα_n1[1:model.ϵSize]
        plasticVars.α = ϵᵖα_n1[model.ϵSize+1:model.ϵSize+model.αSize]
        if algoTangent == true
            ##Calculation of Algorthimic Tangent Tensor
            plasticVars.C .= model.ℂ(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
            plasticVars.D .= model.𝔻(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)
            A = calculateMatrix_A(Δλ, plasticVars, model, params)
            fA .= [model.∂𝒇_∂𝛔(plasticVars.σ_mandel, plasticVars.q,  plasticVars, params);
            model.∂𝒇_∂𝐪(plasticVars.σ_mandel, plasticVars.q,  plasticVars, params)]'*A
            #=Isym = [1.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0
            0.0  0.0  0.0  0.5  0.0  0.0
            0.0  0.0  0.0  0.0  0.5  0.0
            0.0  0.0  0.0  0.0  0.0  0.5]
            Isym = [1.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0
            0.0  0.0  0.0  1.0  0.0  0.0
            0.0  0.0  0.0  0.0  1.0  0.0
            0.0  0.0  0.0  0.0  0.0  1.0]=#
            Isym = getOrder4SymIdentity()
            Θh = [model.𝚯(plasticVars.σ_mandel, plasticVars.q,  plasticVars, params);
            model.𝐡(plasticVars.σ_mandel, plasticVars.q, plasticVars, params)]

            𝐈::Array{Float64, 2}  = [Isym zeros(model.ϵSize, model.αSize); zeros(model.αSize, model.ϵSize) 0.0]
            CTemp::Array{Float64, 2} = A*𝐈 .- A*Θh*(fA*𝐈/(fA*Θh))
            plasticVars.Cᵀ .= get_Pᵀ()*CTemp[1:model.ϵSize, 1:model.ϵSize]*get_P()
        end
    else
        if algoTangent == true
            plasticVars.Cᵀ .= get_Pᵀ()*plasticVars.C*get_P()
        end
    end
    plasticVars.σ_voigt .= get_Pᵀ()*plasticVars.σ_mandel
    return plasticVars.σ_voigt
end

"""This function finds the tangent tensor numerically.

    Cᵀ =  findNumerical_Cᵀ(plasticVars, model, stateDict, params, elementNo, IntegrationPt)
"""

function findNumerical_Cᵀ(plasticVars::PlasticVars, model::PlasticModel,
    params::Parameters, stateDict::Dict{T},
    elementNo::Int64, IntegrationPt::Int64) where T
    stateDictBuffer = createStateDict()
    plasticVarsNew = SmallStrainPlastic.initPlasticVars(model)
    plasticVarsNew.C = plasticVars.C
    Cᵀ = zeros(model.ϵSize,model.ϵSize)
    h = 1e-8
    plasticVarsNew.ϵ = deepcopy(plasticVars.ϵ)
    for i ∈ 1:model.ϵSize
        SmallStrainPlastic.checkPlasticState!(plasticVarsNew, model,
        params, stateDict, stateDictBuffer, elementNo, IntegrationPt)
        σ_old = deepcopy(plasticVarsNew.σ_mandel)
        #println("plasticVarsNew.ϵ = ", plasticVarsNew.ϵ)
        #h = ϵ[i] == 0.0 ? sqrt(eps(1.0)) : sqrt(eps(ϵ[i]))*ϵ[i]
        plasticVarsNew.ϵ[i] +=h
        #println("plasticVarsNew.ϵ + h = ", plasticVarsNew.ϵ)
        SmallStrainPlastic.checkPlasticState!(plasticVarsNew, model,
            params, stateDict, stateDictBuffer, elementNo, IntegrationPt)
        σ_new = deepcopy(plasticVarsNew.σ_mandel)
        Cᵀ[:,i] = (σ_new-σ_old)/h
        #println((σ_new-σ_old))
        plasticVarsNew = SmallStrainPlastic.initPlasticVars(model)
        plasticVarsNew.C = plasticVars.C
        plasticVarsNew.ϵ = deepcopy(plasticVars.ϵ)
    end
    return get_Pᵀ()*Cᵀ*get_P()
end

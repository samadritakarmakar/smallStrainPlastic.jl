
function 𝒇_j2(σ_mandel::Array{T1,1}, q::Array{T2,1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    σ_y::Float64 = params.f
    σ_mandel -= 1.0/3.0*trace(σ_mandel)*getOrder2Identity()
    f = sqrt(3.0/2.0)*frobeniusNorm_p2(σ_mandel)-(σ_y-q[1])
    #println(f)
    return f #<= 0.0 ? 0.0 : f
end

function ∂𝒇_∂𝛔_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    σ_y::Float64 = params.∂f_∂σ
    f(σ) = 𝒇_j2(σ, q, plasticVars, params)
    ∂f_∂σ = ForwardDiff.gradient(f, σ_mandel)
    return ∂f_∂σ
end

function ∂𝒇_∂𝐪_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    ∂f_∂q  = zeros(1)
    ∂f_∂q[1,1] = params.H != 0.0 ? 1.0 : 0.0
    return ∂f_∂q
end

𝚯_j2 = ∂𝒇_∂𝛔_j2

function ∂𝚯_∂𝛔_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    f(σ) = 𝒇_j2(σ, q, plasticVars, params)
    ∂Θ_∂σ = ForwardDiff.hessian(f, σ_mandel)
    return ∂Θ_∂σ
end

function ∂𝚯_∂𝐪_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    ∂Θ_∂q = zeros(9,1)
    return ∂Θ_∂q
end

function 𝐡_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    h = zeros(1)
    h[1] = params.H != 0.0 ? 1.0 : 0.0
    return h
end

#𝐡_j2! = ∂𝒇_∂𝐪_j2!

function ∂𝐡_∂𝛔_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    ∂h_∂σ = zeros(1,9)
    ∂h_∂σ[1,1] = 0.0
    return ∂h_∂σ
end

function ∂𝐡_∂𝐪_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    ∂h_∂q = zeros(1,1)
    ∂h_∂q[1,1] = 0.0
    return ∂h_∂q
end

function 𝓗_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
     α::Array{Float64, 1}, plasticVars::PlasticVars,
    params::ModelParams) where {T1, T2}

    H = zeros(1)
    H[1] = params.H*α[1]
    return H
end

function ℂ_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    ##Do nothing
    return plasticVars.C
end

function 𝔻_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    D = zeros(1,1)
    D[1,1] = params.H != 0.0 ? params.H : 1.0
    return D
end

j2Model = PlasticModel(𝒇_j2, ∂𝒇_∂𝛔_j2, ∂𝒇_∂𝐪_j2, 𝚯_j2,
∂𝚯_∂𝛔_j2, ∂𝚯_∂𝐪_j2, 𝐡_j2, ∂𝐡_∂𝛔_j2, ∂𝐡_∂𝐪_j2,
𝓗_j2, ℂ_j2, 𝔻_j2, 9, 6, 1)

function initParams_j2(σ_y::Float64,  params_H::Float64)
    params_f::Float64 = σ_y
    params_∂f_∂σ::Float64 = σ_y

    return ModelParams{Float64, Float64, Int64, Int64,
        Int64, Int64, Int64, Int64, Int64,
        Float64, Int64, Int64}(params_f, params_∂f_∂σ, 0,0,0,0,
        0, 0,0, params_H,0,0)
end

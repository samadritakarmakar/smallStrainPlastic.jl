
function 𝒇_j2(σ_mandel::Array{T1,1}, q::Array{T2,1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    σ_y::Float64 = params.f
    σ_mandel -= 1.0/3.0*trace(σ_mandel)*getOrder2Identity()
    f = sqrt(3.0/2.0)*frobeniusNorm_p2(σ_mandel-q[2:end])-(σ_y-q[1])
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

    f(q) = 𝒇_j2(σ_mandel, q, plasticVars, params)
    #∂f_∂q[1,1] = params.H != 0.0 ? 1.0 : 0.0
    ∂f_∂q = ForwardDiff.gradient(f, q)
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

    dfdσ(q) =∂𝒇_∂𝛔_j2(σ_mandel, q, plasticVars, params)
    ∂Θ_∂q = ForwardDiff.jacobian(dfdσ, q)
    return ∂Θ_∂q
end

#=function 𝐡_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    h = zeros(1)
    h[1] = params.H != 0.0 ? 1.0 : 0.0
    return h
end=#

𝐡_j2 = ∂𝒇_∂𝐪_j2

function ∂𝐡_∂𝐪_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    f(q) = 𝒇_j2(σ_mandel, q, plasticVars, params)
    ∂h_∂q = ForwardDiff.hessian(f, q)
    return ∂h_∂q
end

#∂𝐡_∂𝛔_j2 = ∂𝚯_∂𝐪_j2

function ∂𝐡_∂𝛔_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    return ∂𝚯_∂𝐪_j2(σ_mandel, q, plasticVars, params)'
end

function ℂ_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    ##Do nothing
    return plasticVars.C
end

function 𝔻_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
    plasticVars::PlasticVars, params::ModelParams) where {T1, T2}

    D = zeros(10,10)
    D[1,1] = params.H[1] != 0.0 ? params.H[1] : 1.0
    D[2:end, 2:end] = params.H[2] != 0.0 ? params.H[2]*plasticVars.C : Array{Float64, 2}(I, 9, 9)
    return D
end

function 𝓗_j2(σ_mandel::Array{T1, 1}, q::Array{T2, 1},
     α::Array{Float64, 1}, plasticVars::PlasticVars,
    params::ModelParams) where {T1, T2}

    H = 𝔻_j2(σ_mandel, q, plasticVars, params)*α
    return H
end

j2Model = PlasticModel(𝒇_j2, ∂𝒇_∂𝛔_j2, ∂𝒇_∂𝐪_j2, 𝚯_j2,
∂𝚯_∂𝛔_j2, ∂𝚯_∂𝐪_j2, 𝐡_j2, ∂𝐡_∂𝛔_j2, ∂𝐡_∂𝐪_j2,
𝓗_j2, ℂ_j2, 𝔻_j2, 9, 6, 10)

function initParams_j2(σ_y::Float64,  params_Hi::Float64, kinematicFraction::Float64 = 0.0)
    params_f::Float64 = σ_y
    params_∂f_∂σ::Float64 = σ_y
    params_H::Tuple{Float64, Float64} = (params_Hi, kinematicFraction)

    return ModelParams{Float64, Float64, Int64, Int64,
        Int64, Int64, Int64, Int64, Int64,
        Tuple{Float64, Float64}, Int64, Int64}(params_f, params_∂f_∂σ, 0,0,0,0,
        0, 0,0, params_H, 0,0)
end

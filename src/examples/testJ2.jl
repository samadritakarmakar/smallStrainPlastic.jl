using SmallStrainPlastic, LinearAlgebra, Tensors, Plots

function testJ2()
    σ_y = 200.0
    E = 200.0e3
    ν = 0.3
    stateDict = createStateDict()
    stateDictBuffer = createStateDict()
    plasticVars = SmallStrainPlastic.initPlasticVars(SmallStrainPlastic.j2Model)
    #plasticVars.C = SmallStrainPlastic.createVoigtElasticTensor(E, ν)
    plasticVars.C .= SmallStrainPlastic.getMandelElasticTensor(E, ν)

    params_J2 = SmallStrainPlastic.initParams_j2(σ_y, 20.0e3)

    σₘArray::Array{Float64, 1} = zeros(0)
    𝐬Array::Array{Float64, 1} = zeros(0)
    ϵₘArray::Array{Float64, 1} = zeros(0)
    𝒆Array::Array{Float64, 1} = zeros(0)
    iArray::Array{Int64, 1} = zeros(0)
    #plasticVars.ϵ[1] += 20e-4
    for i ∈ 1:20
        if (i<=20)
            plasticVars.ϵ[1] += 1e-4
        elseif (i>20 && i<=55)
            plasticVars.ϵ[1] -= 1e-4
        else
            plasticVars.ϵ[1] += 1e-4
        end
        SmallStrainPlastic.checkPlasticState!(plasticVars, SmallStrainPlastic.j2Model,
        params_J2, stateDict, stateDictBuffer, 1, 1, algoTangent = true)
        println("Cᵀ Algorithmic", plasticVars.Cᵀ)
        #println(" ϵᵖ = ", plasticVars.ϵᵖ, " α = ", plasticVars.α, " σ = ", plasticVars.σ_voigt)
        Cᵀ = SmallStrainPlastic.findNumerical_Cᵀ(plasticVars, SmallStrainPlastic.j2Model, params_J2, stateDict, 1, 1)
        #=testϵᵖ = zeros(6)
        #testα = zeros(1)
        #SmallStrainPlastic.getState!(testϵᵖ, testα, stateDict, 1,1)
        #println("testϵᵖ = ", testϵᵖ, "\ntestα = ", testα)=#
        println("Cᵀ Numerical", Cᵀ)
        σₘ, 𝐬 = SmallStrainPlastic.get_σₘ_𝐬_mandel(plasticVars.σ_mandel)
        push!(σₘArray, σₘ)
        push!(𝐬Array, 𝐬)
        ϵₘ, 𝒆 = get_ϵₘ_𝒆_mandel(plasticVars.ϵ)
        push!(ϵₘArray, ϵₘ)
        push!(𝒆Array, 𝒆)
        #push!(iArray, i)
        SmallStrainPlastic.updateStateDict4rmBuffer!(stateDict, stateDictBuffer)
    end

    plot(ϵₘArray, 𝐬Array, legend=false)#, seriestype = :scatter)
end

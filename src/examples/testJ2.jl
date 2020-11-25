using SmallStrainPlastic, Plots

function testJ2()
    σ_y = 200.0
    E = 200e3
    ν = 0.3
    plasticVars = SmallStrainPlastic.initPlasticVars(SmallStrainPlastic.j2Model)
    plasticVars.C = SmallStrainPlastic.createVoigtElasticTensor(E, ν)
    params_J2 = SmallStrainPlastic.initParams_j2(σ_y, -20.0e3)
    σₘArray::Array{Float64, 1} = zeros(0)
    𝐬Array::Array{Float64, 1} = zeros(0)
    ϵₘArray::Array{Float64, 1} = zeros(0)
    𝒆Array::Array{Float64, 1} = zeros(0)
    iArray::Array{Int64, 1} = zeros(0)
    SmallStrainPlastic.tolerance.R = 1e-8
    SmallStrainPlastic.tolerance.maxIter = 2
    for i ∈ 1:90
        if (i<=20)
            plasticVars.ϵ[1] += 1e-4
        elseif (i>20 && i<=55)
            plasticVars.ϵ[1] -= 1e-4
        else
            plasticVars.ϵ[1] += 1e-4
        end
        SmallStrainPlastic.checkPlasticState!(plasticVars, SmallStrainPlastic.j2Model, params_J2, 1, 1)
        println(" ϵᵖ = ", plasticVars.ϵᵖ, " α = ", plasticVars.α)
        σₘ, 𝐬 = SmallStrainPlastic.get_σₘ_𝐬(plasticVars.σ_voigt)
        push!(σₘArray, σₘ)
        push!(𝐬Array, 𝐬)
        ϵₘ, 𝒆 = get_ϵₘ_𝒆(plasticVars.ϵ)
        push!(ϵₘArray, ϵₘ)
        push!(𝒆Array, 𝒆)
        push!(iArray, i)
        SmallStrainPlastic.updateStateDict4rmBuffer()
    end

    plot(ϵₘArray, 𝐬Array, legend=false)#, seriestype = :scatter)
end

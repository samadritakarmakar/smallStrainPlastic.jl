using SmallStrainPlastic, Plots

function testJ2()
    σ_y = 200.0
    E = 200e3
    ν = 0.3
    plasticVars = SmallStrainPlastic.initPlasticVars(SmallStrainPlastic.j2Model)
    plasticVars.C = SmallStrainPlastic.createVoigtElasticTensor(E, ν)
    params_J2 = SmallStrainPlastic.initParams_j2(σ_y, 0.0)
    𝒑Array::Array{Float64, 1} = zeros(0)
    𝒒Array::Array{Float64, 1} = zeros(0)
    𝒆Array::Array{Float64, 1} = zeros(0)
    𝒆ₛArray::Array{Float64, 1} = zeros(0)
  iArray::Array{Int64, 1} = zeros(0)
    for i ∈ 1:82

        println("ϵ = ", plasticVars.ϵ, " ϵᵖ = ", plasticVars.ϵᵖ, " α = ", plasticVars.α)
        if (i<=20)
            plasticVars.ϵ[1] += 1e-4
        elseif (i>20 && i<=55)
            println("runs!")
            plasticVars.ϵ[1] -= 1e-4
        else
            plasticVars.ϵ[1] += 1e-4
        end
        SmallStrainPlastic.checkPlasticState!(plasticVars, SmallStrainPlastic.j2Model, params_J2, 1, 1)
        𝒑, 𝒒 = SmallStrainPlastic.get_𝒑_𝒒(plasticVars.σ_voigt)
        push!(𝒑Array, 𝒑)
        push!(𝒒Array, 𝒒)
        𝒆, 𝒆ₛ = get_𝒆_𝒆ₛ(plasticVars.ϵ)
        push!(𝒆Array, 𝒆)
        push!(𝒆ₛArray, 𝒆ₛ)
        push!(iArray, i)
    end
    plot(𝒆Array, 𝒒Array, legend=false)#, seriestype = :scatter)
end

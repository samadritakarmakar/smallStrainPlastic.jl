using SmallStrainPlastic, LinearAlgebra, Tensors, Plots

function testJ2()
    σ_y = 200.0
    E = 200e3
    ν = 0.3
    plasticVars = SmallStrainPlastic.initPlasticVars(SmallStrainPlastic.j2Model)
    plasticVars.C = SmallStrainPlastic.createVoigtElasticTensor(E, ν)

    params_J2 = SmallStrainPlastic.initParams_j2(σ_y, 20e3)

    σₘArray::Array{Float64, 1} = zeros(0)
    𝐬Array::Array{Float64, 1} = zeros(0)
    ϵₘArray::Array{Float64, 1} = zeros(0)
    𝒆Array::Array{Float64, 1} = zeros(0)
    iArray::Array{Int64, 1} = zeros(0)
    #σ_buffer = zeros(6)
    #ϵ_buffer = zeros(6)
    for i ∈ 1:820
        if (i<=200)
            plasticVars.ϵ[1] += 1e-5
        elseif (i>200 && i<=550)
            plasticVars.ϵ[1] -= 1e-5
        else
            plasticVars.ϵ[1] += 1e-5
        end
        SmallStrainPlastic.checkPlasticState!(plasticVars, SmallStrainPlastic.j2Model, params_J2, 1, 1)
        #println(" ϵᵖ = ", plasticVars.ϵᵖ, " α = ", plasticVars.α)
        Cᵀ::SymmetricTensor{4,3} = Tensors.fromvoigt(SymmetricTensor{4,3},plasticVars.Cᵀ)
        #println("norm(Cᵀ) = ", norm(Cᵀ), " rank(plasticVars.Cᵀ) ", rank(plasticVars.Cᵀ))
        #Δσ = plasticVars.σ_voigt - σ_buffer
        #Δϵ = plasticVars.ϵ - ϵ_buffer
        #σ_buffer = deepcopy(plasticVars.σ_voigt)
        #ϵ_buffer = deepcopy(plasticVars.ϵ)
        #println(Δσ./Δϵ')
        #println(plasticVars.Cᵀ)
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

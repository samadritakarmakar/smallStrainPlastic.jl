using SmallStrainPlastic

function testJ2()
    σ_y = 200.0
    E = 200e3
    ν = 0.3
    plasticVars = SmallStrainPlastic.initPlasticVars(SmallStrainPlastic.j2Model)
    plasticVars.C = SmallStrainPlastic.createVoigtElasticTensor(E, ν)
    params_J2 = SmallStrainPlastic.initParams_j2(σ_y, 0.0)
    for i ∈ 1:20
        SmallStrainPlastic.checkPlasticState!(plasticVars, SmallStrainPlastic.j2Model, params_J2, 1, 1)
        println("ϵ = ", plasticVars.ϵ, " ϵᵖ = ", plasticVars.ϵᵖ, " σ_voigt = ", plasticVars.σ_voigt)
        plasticVars.ϵ[1] += 1e-3
    end
    #println("C = ", C)
end

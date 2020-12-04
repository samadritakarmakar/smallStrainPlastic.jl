module SmallStrainPlastic
using LinearAlgebra, Tensors
include("plasticity/plasticModel.jl")
include("plasticity/returnMapping.jl")
include("plasticity/state.jl")
include("tools/denseJacobian.jl")
include("tools/TensorTools.jl")
include("plasticModels/j2.jl")



#From plasticity/plasticModel.jl
export PlasticModel, PlasticVars, initPlasticVars
export Parameters, ModelParams, tolerance
#tolerance = Tolerance(1e-8, 1e-8, 1000)
#From plasticity
##from state.jl
export State, getState!, updateStateDict!, createStateDict, updateStateDict4rmBuffer!
#stateDict = createStateDict()
#stateDictBuffer = createStateDict()
##from returnMapping.jl
export checkPlasticState!

#From tools/
##from denseJacobian
export denseJacobian!, denseJacobian
##from TensorTools
export createVoigtElasticTensor, getProjectionTensor4, get_σₘ_𝐬, get_ϵₘ_𝒆

#From plasticModels/
##from j2.jl
export 𝒇_j2, ∂𝒇_∂𝛔_j2!, ∂𝒇_∂𝐪_j2!, ∂Θ_∂𝛔_j2!, ∂Θ_∂𝐪_j2!
export ∂𝐡_∂𝛔_j2!, ∂𝐡_∂𝐪_j2!, 𝓗_j2!, ℂ_j2!, 𝔻_j2!
export j2Model, initParams_j2
end # module

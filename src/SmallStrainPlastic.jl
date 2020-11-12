module SmallStrainPlastic
using LinearAlgebra
include("plasticity/plasticModel.jl")
include("plasticity/returnMapping.jl")
include("plasticity/state.jl")
include("tools/denseJacobian.jl")

#From plasticity/plasticModel.jl
export PlasticModel

#From plasticity/state.jl
export State, getState!, updateStateDict!, createStateDict

#From tools/denseJacobian
export denseJacobian!, denseJacobian
end # module

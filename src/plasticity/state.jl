"""This structure saves the state of the material, meaning it's plastic strain and hardening variable."""
mutable struct State
    ϵᵖ::Array{Float64,1}
    α::Array{Float64,1}
end

"""
This function gets the state of the material, meaning it's plastic strain and hardening variable.
If they exist in the Dictionary for the given material/integration point in the given element,
it updates the data with the available data in stateDict.
If they don't exist, it just fills the state varibles with zeros.

    getState!(ϵᵖ, α, stateDict, elementNo, integrationPt)
"""
function getState!(ϵᵖ::Array{Float64,1}, α::Array{Float64,1}, stateDict::Dict{Tuple{Int64, Int64}, State},
    elementNo::Int64= 1, integrationPt::Int64=1)

    if (elementNo, integrationPt) ∈ keys(stateDict)
        ϵᵖ = stateDict[elementNo, integrationPt].ϵᵖ
        α = stateDict[elementNo, integrationPt].α
    else
        fill!(ϵᵖ, 0.0)
        fill!(α, 0.0)
    end
    return nothing
end

"""
This function updates the StateDict according to the passed data of ϵᵖ and α for a specific element number and
an integration point within the given element.

    updateStateDict!(ϵᵖ, α, stateDict, elementNo, integrationPt)
"""
function updateStateDict!(ϵᵖ::Array{Float64,1}, α::Array{Float64,1}, stateDict::Dict{Tuple{Int64, Int64}, State},
    elementNo::Int64= 1, integrationPt::Int64=1)
    stateDict[elementNo, integrationPt] = State(ϵᵖ, α)
    return nothing
end

"""
This function creates a Dictionary of type State to store the state of the material, meaning it's plastic strain and hardening variable.

    stateDict = createStateDict()
"""
function createStateDict()
    return Dict{Tuple{Int64, Int64}, State}()
end

function updateStateDict!(stateDict::Dict{Tuple{Int64, Int64}, State},
    stateDictBuffer::Dict{Tuple{Int64, Int64}, State})
    buffKeys = keys(SmallStrainPlastic.stateDictBuffer)
    #println("buffKeys = ", buffKeys)
    for buffKey ∈ collect(buffKeys)
        stateDict[buffKey...] = stateDictBuffer[buffKey...]
    end
end

function updateStateDict4rmBuffer()
    updateStateDict!(stateDict, stateDictBuffer)
end

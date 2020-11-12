"""This structure saves the state of the material, meaning it's plastic strain and hardening variable."""
mutable struct State
    𝛆ᵖ::Array{Float64,1}
    𝛂::Array{Float64,1}
end

"""
This function gets the state of the material, meaning it's plastic strain and hardening variable.
If they exist in the Dictionary for the given material/integration point in the given element,
it updates the data with the available data in stateDict.
If they don't exist, it just fills the state varibles with zeros.

    getState!(𝛆ᵖ, 𝛂, stateDict, elementNo, integrationPt)
"""
function getState!(𝛆ᵖ::Array{Float64,1}, 𝛂::Array{Float64,1}, stateDict::Dict{Tuple{Int64, Int64}, State},
    elementNo::Int64= 1, integrationPt::Int64=1)
    
    if (elementNo, integrationPt) ∈ keys(stateDict)
        𝛆ᵖ = stateDict[elementNo, integrationPt].𝛆ᵖ
        𝛂 = stateDict[elementNo, integrationPt].𝛂
    else
        fill!(materialState.𝛆ᵖ, 0.0)
        fill!(materialState.𝛂, 0.0)
    end
    return nothing
end

"""
This function updates the StateDict according to the passed data of 𝛆ᵖ and 𝛂 for a specific element number and
an integration point within the given element.

    updateStateDict!(𝛆ᵖ, 𝛂, stateDict, elementNo, integrationPtNo)
"""
function updateStateDict!(𝛆ᵖ::Array{Float64,1}, 𝛂::Array{Float64,1}, stateDict::Dict{Tuple{Int64, Int64}, State},
    elementNo::Int64= 1, integrationPtNo::Int64=1)
    stateDict[elementNo, integrationPt] = State(𝛆ᵖ, 𝛂)
    return nothing
end

"""
This function creates a Dictionary of type State to store the state of the material, meaning it's plastic strain and hardening variable.

    stateDict = createStateDict()
"""
function createStateDict()
    return Dict{Tuple{Int64, Int64}, State}()
end

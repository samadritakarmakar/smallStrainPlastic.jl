"""A plasticity new plasticity model can be defined by defininf the following functions:

Yield Function:- 𝒇!

The partial of Yield Function with respect to stress, ∂𝒇/∂𝛔:- ∂𝒇_∂𝛔!

The partial of Yield Function with respect to hardening variable, ∂𝒇/∂𝐪: ∂𝒇_∂𝐪!

The function in the flow rule for the plastic strain, 𝛆̇ᵖ = λ̇  Θ(𝛔, 𝐪) :- Θ! #Removed

The partial of plastic strain flow rule function with respect to stress, ∂Θ/∂𝛔:- ∂Θ_∂𝛔!

The partial of plastic strain flow rule function with respect to hardening, ∂Θ/∂𝐪:- ∂Θ_∂𝐪!

The function in the flow rule for the internal variable 𝛂̇ = λ̇  𝐡(𝛔, 𝐪):- 𝐡! #Removed

The partial of plastic strain flow rule function with respect to stress, ∂𝐡/∂𝛔:- ∂𝐡_∂𝛔!

The partial of plastic strain flow rule function with respect to hardening, ∂𝐡/∂𝐪:- ∂𝐡_∂𝐪!

If the evolution of the hardening variable 𝐪̇ is defined as 𝐪̇ = -𝓗(𝛂), then the function it is
dependent on can be written as:- 𝓗! #Removed

If the hardening variable 𝐪̇ is defined as 𝐪̇ = -𝓗(𝛂), then an equivalent to stiffness tensor
ℂ = ∂𝛔/∂𝛆ᵉ can defined as 𝔻 = -∂𝐪/∂𝛂 = ∂𝓗(𝛂)/∂𝐪 :- 𝔻!

For ease of use, defining a function that saves the stiffness tensor is also made available :- ℂ
 ℂ!
"""
struct PlasticModel
    𝒇::Function
    ∂𝒇_∂𝛔!::Function
    ∂𝒇_∂𝐪!::Function
    #Θ!::Function ##Does not seem to be needed here for return mapping algorithm
    ∂Θ_∂𝛔!::Function
    ∂Θ_∂𝐪!::Function
    #𝐡!::Function ##Does not seem to be needed here for return mapping algorithm
    ∂𝐡_∂𝛔!::Function
    ∂𝐡_∂𝐪!::Function
    #𝓗!::Function ##Does not seem to be needed here for return mapping algorithm
    𝔻!::Function
    ℂ!::Function
end

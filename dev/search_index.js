var documenterSearchIndex = {"docs":
[{"location":"#SmallStrainPlastic-Documentation","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"","category":"section"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"SmallStrainPlastic is a library that aims to simplify the addition of different plastic models. The usage may be for a material point to judge the behaviour of a plastic model or to in a finite element model to simulate the behaviour of a plastic material.","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"The assumptions made while developing the library are the following:","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"The flow rule for plastic strain is considered to be:","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"𝛆̇ᵖ = λ̇  Θ(𝛔, 𝐪)","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"The flow rule for internal variable is considered to be:","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"𝛂̇ = λ̇  𝐡(𝛔, 𝐪)","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"The stress is considered as:","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"𝛔 = ℂ:(𝛆 - 𝛆ᵖ)","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"The hardening variable is considered as:","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"𝐪 = - 𝓗(𝛂)","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"Special Unicode characters using eg: \"\\bb\", \"\\bf\", or \"\\bsrc*\" are used to define functions. It is highly recommended such unicode characters are avoided when defining internal variables to avoid confusion.","category":"page"},{"location":"#Plasticity-Model","page":"SmallStrainPlastic Documentation","title":"Plasticity Model","text":"","category":"section"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"\tSmallStrainPlastic.PlasticModel","category":"page"},{"location":"#SmallStrainPlastic.PlasticModel","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic.PlasticModel","text":"A plasticity new plasticity model can be defined by defining the following functions:\n\nYield Function:-\n\n𝒇(σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)\n\nThe partial of Yield Function with respect to stress, ∂𝒇/∂𝛔:-\n\n∂𝒇_∂𝛔!(∂f_∂σ::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)\n\nThe partial of Yield Function with respect to hardening variable, ∂𝒇/∂𝐪:\n\n∂𝒇_∂𝐪!(∂f_∂q::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1},  params::ModelParams)\n\nThe function in the flow rule for the plastic strain, 𝛆̇ᵖ = λ̇  𝚯(𝛔, 𝐪) :-\n\n𝚯!(Θ::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)\n\nThe partial of plastic strain flow rule function with respect to stress, ∂𝚯/∂𝛔:-\n\n∂𝚯_∂𝛔!(∂Θ_∂σ::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)\n\nThe partial of plastic strain flow rule function with respect to hardening, ∂𝚯/∂𝐪:-\n\n∂𝚯_∂𝐪!(∂Θ_∂q::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)\n\nThe function in the flow rule for the internal variable 𝛂̇ = λ̇  𝐡(𝛔, 𝐪):-\n\n𝐡!(h::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)\n\nThe partial of plastic strain flow rule function with respect to stress, ∂𝐡/∂𝛔:-\n\n∂𝐡_∂𝛔!(∂h_∂σ::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, params::ModelParams)\n\nThe partial of plastic strain flow rule function with respect to hardening, ∂𝐡/∂𝐪:-\n\n∂𝐡_∂𝐪!(∂h_∂q::Array{Float64, 2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1},  params::ModelParams)\n\nIf the evolution of the hardening variable 𝐪̇ is defined as 𝐪̇ = -𝓗(𝛂), then the function it is dependent on can be written as:-\n\n𝓗!(H::Array{Float64, 1}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1}, α::Array{Float64, 1}, params::ModelParams)\n\nFor ease of use, defining a function that saves the stiffness tensor is also made available :-\n\nℂ!(C::Array{Float64,2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1},  params::ModelParams)\n\nIf the hardening variable 𝐪̇ is defined as 𝐪̇ = -𝓗(𝛂), then an equivalent to stiffness tensor defined as ℂ = ∂𝛔/∂𝛆ᵉ, we can defined as 𝔻 = -∂𝐪/∂𝛂 = ∂𝓗(𝛂)/∂𝐪 :-\n\n𝔻!(D::Array{Float64,2}, σ_voigt::Array{Float64, 1}, q::Array{Float64, 1},  params::ModelParams)\n\n\n\n\n\n","category":"type"},{"location":"#State-of-Plasticity","page":"SmallStrainPlastic Documentation","title":"State of Plasticity","text":"","category":"section"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"\tSmallStrainPlastic.State","category":"page"},{"location":"#SmallStrainPlastic.State","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic.State","text":"This structure saves the state of the material, meaning it's plastic strain and hardening variable.\n\n\n\n\n\n","category":"type"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"\tSmallStrainPlastic.createStateDict","category":"page"},{"location":"#SmallStrainPlastic.createStateDict","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic.createStateDict","text":"This function creates a Dictionary of type State to store the state of the material, meaning it's plastic strain and hardening variable.\n\nstateDict = createStateDict()\n\n\n\n\n\n","category":"function"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"\tSmallStrainPlastic.updateStateDict!","category":"page"},{"location":"#SmallStrainPlastic.updateStateDict!","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic.updateStateDict!","text":"This function updates the StateDict according to the passed data of ϵᵖ and α for a specific element number and an integration point within the given element.\n\nupdateStateDict!(ϵᵖ, α, stateDict, elementNo, integrationPt)\n\n\n\n\n\n","category":"function"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"\tSmallStrainPlastic.getState!","category":"page"},{"location":"#SmallStrainPlastic.getState!","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic.getState!","text":"This function gets the state of the material, meaning it's plastic strain and hardening variable. If they exist in the Dictionary for the given material/integration point in the given element, it updates the data with the available data in stateDict. If they don't exist, it just fills the state varibles with zeros.\n\ngetState!(ϵᵖ, α, stateDict, elementNo, integrationPt)\n\n\n\n\n\n","category":"function"},{"location":"#Easier-Finding-of-Jacobians","page":"SmallStrainPlastic Documentation","title":"Easier Finding of Jacobians","text":"","category":"section"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"\tSmallStrainPlastic.denseJacobian!","category":"page"},{"location":"#SmallStrainPlastic.denseJacobian!","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic.denseJacobian!","text":"When working with plastic models it many times it becomes difficult to find the analytically jacobian of a vector function. With denseJacobian we try to make that easy to do. It uses finite difference to do this. Finite Differences are far from the best solution to find jacobians but are the easiest to implement. Hence this solution.\n\ndenseJacobian!(jacobian::Array{Float64,2}, f::Function, x::Array{Float64,1})\n\nHere \"jacobian\" must have row size equal to length of vector returned by Function \"f\" and column size equal to length of vector \"x\"\n\n\n\n\n\n","category":"function"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"\tSmallStrainPlastic.denseJacobian","category":"page"},{"location":"#SmallStrainPlastic.denseJacobian","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic.denseJacobian","text":"denseJacobian is an easier to use version of denseJacobian!. This is less efficient than denseJacobian!. Hence whenever possible denseJacobian! must be used.\n\n\n\n\n\n","category":"function"},{"location":"#Return-Mapping-Alogrithm","page":"SmallStrainPlastic Documentation","title":"Return Mapping Alogrithm","text":"","category":"section"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"\tSmallStrainPlastic.returnMapping!","category":"page"},{"location":"#SmallStrainPlastic.returnMapping!","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic.returnMapping!","text":"This function is responsible for executing the return mapping algorithm. It does so by calculating the evolution of the plastic strain using the Closest Point Projection method. The following formulations are used\n\nHere's an equation:\n\nd(Deltalambda) = fracf^k - beginbmatrixpartial f^kpartial sigma  partial f^k partial qendbmatrixbeginbmatrixAendbmatrixbeginbmatrix R endbmatrix beginbmatrixpartial f^kpartial sigma  partial f^k partial qendbmatrixbeginbmatrixAendbmatrixbeginbmatrixTheta  h endbmatrix\n\nwhere:\n\nbeginbmatrixRendbmatrix = -beginbmatrix epsilon^p_n+1  alpha_n+1 endbmatrix +beginbmatrix epsilon^p_n  alpha_n endbmatrix +Deltalambdabeginbmatrix Theta(sigma_n+1 q_n+1)  h(sigma_n+1 q_n+1) endbmatrix\n\nbeginbmatrixAendbmatrix^-1 = beginbmatrix bmC^-1 + Deltalambda fracpartial Thetapartialsigma_n+1  Deltalambda fracpartial Thetapartial q_n+1  Deltalambda fracpartial hpartialsigma_n+1  bmD^-1 + Deltalambda fracpartial hpartial q_n+1 endbmatrix\n\nThe Strain 𝛆ᵖ and the internal variable 𝛂 are updated as,\n\nbeginbmatrixDelta epsilon^p  Delta alpha endbmatrix = beginbmatrixbmC^-1  0  0  bmD^-1 endbmatrix beginbmatrixAendbmatrix beginbmatrixTheta  h endbmatrix d(Deltalambda)\n\n\n\n\n\n","category":"function"},{"location":"#An-example-using-J2-Plastic-Model","page":"SmallStrainPlastic Documentation","title":"An example using J2 Plastic Model","text":"","category":"section"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"using SmallStrainPlastic, Plots\nfunction testJ2()\n\tσ_y = 200.0\n\tE = 200e3\n\tν = 0.3\n\tplasticVars =SmallStrainPlastic.initPlasticVars(SmallStrainPlastic.j2Model)\n  \tplasticVars.C = SmallStrainPlastic.createVoigtElasticTensor(E, ν)\n\tparams_J2 = SmallStrainPlastic.initParams_j2(σ_y, 0.0)\n\t𝒑Array::Array{Float64, 1} = zeros(0)\n\t𝒒Array::Array{Float64, 1} = zeros(0)\n\t𝒆Array::Array{Float64, 1} = zeros(0)\n\t𝒆ₛArray::Array{Float64, 1} = zeros(0)\n\tfor i ∈ 1:82\n    \t\tif (i<=20)\n        \t\tplasticVars.ϵ[1] += 1e-4\n    \t\telseif (i>20 && i<=55)\n        \t\tplasticVars.ϵ[1] -= 1e-4\n    \t\telse\n       \t\t \tplasticVars.ϵ[1] += 1e-4\n    \t\tend\n    \t\tSmallStrainPlastic.checkPlasticState!(plasticVars, SmallStrainPlastic.j2Model, params_J2, 1, 1)\n   \t\t 𝒑, 𝒒 = SmallStrainPlastic.get_𝒑_𝒒(plasticVars.σ_voigt)\n    \t\tpush!(𝒑Array, 𝒑)\n    \t\tpush!(𝒒Array, 𝒒)\n    \t\t𝒆, 𝒆ₛ = get_𝒆_𝒆ₛ(plasticVars.ϵ)\n    \t\tpush!(𝒆Array, 𝒆)\n   \t\tpush!(𝒆ₛArray, 𝒆ₛ) \n\tend\n\tplot(𝒆Array, 𝒒Array, legend=false)#, seriestype = :scatter)\nend","category":"page"},{"location":"","page":"SmallStrainPlastic Documentation","title":"SmallStrainPlastic Documentation","text":"You should get a plot like this: (Image: Plot Perfect Plasticity)","category":"page"}]
}

function createVoigtElasticTensor(E::Float64, ν::Float64)::Array{Float64, 2}
    c::Float64 = E/((1+ν)*(1-2*ν))
    C::Array{Float64, 2} = zeros(6,6)
    C = [1-ν ν ν 0 0 0;
         ν 1-ν ν 0 0 0;
         ν ν 1-ν 0 0 0;
         0 0 0 (1-2*ν)/2 0 0;
         0 0 0 0 (1-2*ν)/2 0;
         0 0 0 0 0 (1-2*ν)/2]
    C = c*C
    return C
end

"""From the book The Finite Element Method for Solid and Structural Mechanics,
Seventh Edition, O.C. Zienkiewicz, R.L. Taylor, D.D. Fox.
Pᵀ = [2.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  2.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  2.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0]
 and,
 σ_voigt = Pᵀ*σ_Mandel
 ϵ_Mandel = P*ϵ_voigt
"""
function get_Pᵀ()
    Pᵀ= zeros(6,9)
    Pᵀ[1:3,1:3]= Array{Float64, 2}(I, 3, 3)
    col::Int64 = 4
    for row ∈ 4:6
        Pᵀ[row,col:col+1] = 0.5*ones(2)
        col +=2
    end
    return Pᵀ
end

"""From the book The Finite Element Method for Solid and Structural Mechanics,
Seventh Edition, O.C. Zienkiewicz, R.L. Taylor, D.D. Fox.
Pᵀ = [2.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  2.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  2.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0]
 and,
 σ_voigt = Pᵀ*σ_Mandel
 ϵ_Mandel = P*ϵ_voigt
"""
function get_P()
    return get_Pᵀ()'
end

function getTensorMapping()::Dict{Int64, Int64}
    mapDict::Dict{Int64, Int64} = Dict{Int64, Int64}()
    mapDict[11] = 1
    mapDict[22] = 2
    mapDict[33] = 3
    mapDict[12] = 4
    mapDict[21] = 5
    mapDict[23] = 6
    mapDict[32] = 7
    mapDict[31] = 8
    mapDict[13] = 9
    return mapDict
end

function getVoigtIndex(mapDict::Dict{Int64, Int64}, i::Int64, j::Int64)::Int64
    return mapDict[10*i+j]
end

"""Double contraction of two tensors using Mandel Notation as per "The Finite Element Method for Solid and Structural Mechanics,
Seventh Edition, O.C. Zienkiewicz, R.L. Taylor, D.D. Fox."

result = doubleContract(Array1, Array2)
"""
function doubleContract(Array1::Array{T,1}, Array2::Array{T,1}) where T
    mapDict = getTensorMapping()
    if (size(Array1, 2) ==1)
        Array1 = Array1' #Transposes Array1 so that the use Array1[ij,mn] is still valid
    end
    rows = size(Array1, 1)
    cols = size(Array2, 2)
    @assert ((rows == 1 || rows == 9) && (cols ==1 || cols ==9)) "Array1 or Array2
    must either have a size of 9x1 or 9x9 representing a second order tensor or a
    4th order tensor respectively."
    ijLength = rows > 1 ? 3 : 1
    klLength = cols > 1 ? 3 : 1
    #R_{ijkl} = C_{ijmn}D_{opkl}
    R = zeros(T, rows, cols)
    for k ∈ 1:klLength
        for l ∈ 1:klLength
            kl = getVoigtIndex(mapDict, k,l)
            for i ∈ 1:ijLength
                for j ∈ 1:ijLength
                    ij = getVoigtIndex(mapDict, i,j)
                    for m ∈ 1:3
                        for n ∈ 1:3
                            mn = getVoigtIndex(mapDict, m,n)
                            R[ij,kl] += Array1[ij,mn]*Array2[mn,kl]
                        end
                    end
                end
            end
        end
    end
    return R
end
"""2nd Order tensor Norm save in array as per  mandel notation given in "The Finite Element Method for Solid and Structural Mechanics,
Seventh Edition, O.C. Zienkiewicz, R.L. Taylor, D.D. Fox."

    frobeniusNorm_p2(array)
"""
function doubleContract(array::Array{T,1}) where T
    fN = doubleContract(array, array)
    return sqrt(fN[1])
end

"""Trace of second order mandel notation tensor.

    arrayTrace = trace(array)
"""
function trace(array::Array{T,1}) where T
    m = zeros(9)
    m[1:3] = ones(3)
    return m'*array
end

"""Generates 2nd Order Identity tensor in mandel notation

    Order2I = getOrder2Identity()
"""
function getOrder2Identity()
    δ(i, j) = i == j ? 1.0 : 0.0
    mapDict = getTensorMapping()
    I = zeros(9)
    for i ∈ 1:3
        for j ∈ 1:3
            ij = getVoigtIndex(mapDict, i, j)
            I[ij] += δ(i,j)
        end
    end
    return I
end


"""Generates 4th Order Identity tensor in mandel notation

    Order4I = getOrder4Identity()
"""
function getOrder4Identity()
    δ(i, j) = i == j ? 1.0 : 0.0
    mapDict = getTensorMapping()
    I = zeros(9,9)
    for k ∈ 1:3
        for l ∈ 1:3
            kl = getVoigtIndex(mapDict, k, l)
            for i ∈ 1:3
                for j ∈ 1:3
                    ij = getVoigtIndex(mapDict, i, j)
                    I[ij,kl] += δ(i,j)*δ(k,l)
                end
            end
        end
    end
    return I
end


"""Generates 4th Order Symmetric Identity tensor in mandel notation

    Order4ISym = getOrder4SymIdentity()
"""
function getOrder4SymIdentity()
    δ(i, j) = i == j ? 1.0 : 0.0
    mapDict = getTensorMapping()
    I = zeros(9,9)
    for k ∈ 1:3
        for l ∈ 1:3
            kl = getVoigtIndex(mapDict, k, l)
            for i ∈ 1:3
                for j ∈ 1:3
                    ij = getVoigtIndex(mapDict, i, j)
                    I[ij,kl] +=0.5*(δ(i,k)*δ(j,l)+δ(i,l)*δ(j,k))
                end
            end
        end
    end
    return I
end

"""Generates Elastic 4th order Tensor as per Mandel Notation mentioned in "The Finite Element Method for Solid and Structural Mechanics,
Seventh Edition, O.C. Zienkiewicz, R.L. Taylor, D.D. Fox."

    C_mandel = getMandelElasticTensor(200e3, 0.3)
"""
function getMandelElasticTensor(E::Float64, ν::Float64)::Array{Float64, 2}
    λ = ν*E/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    I4 = getOrder4Identity()
    #I4Sym = getOrder4SymIdentity()
    #return λ*I4 .+ 2*μ*I4Sym
    return λ*I4 .+ 2*μ*Array{Float64,2}(I, 9,9)
end

"""Transform mandel notation 2nd or 4th order tensor to voigt notation

    arrayVoigt = mandel2voigt(arrayMandel)
"""
function mandel2voigt(array::Array{T,1}) where T
    rowSize = size(array, 1)
    colSize = size(array, 2)
    @assert (rowSize==9 || colSize ==9) "mandel2voigt function can only be used for arrays of 9x1 or 9x9 size."
    rowLength = rowSize > 1 ? 3 : 1
    colLength = colSize > 1 ? 3 : 1
    Pᵀ = get_Pᵀ()
    if (rowLength ==3 && colLength ==3)
        P = get_P()
        arrayVoigt = Pᵀ*array*P
    else
        arrayVoigt = Pᵀ*array
    end
    return arrayVoigt
end


"""Transform mandel notation stress to voigt notation or engineering stress

    arrayVoigt = getVoigtEngineeringStress(arrayMandel)
"""

function getVoigtEngineeringStress(σ_mandel::Array{T,1}) where T
    return Pᵀ*σ_mandel
end

"""Transform voigt notation strain to  mandel continuum strain

    strainEngineering = mandel2voigt(arrayMandel)
"""
function getContinuumMandelStrain(ϵ_MandelContinuum::Array{T,1}) where T
    return get_P()*ϵ_MandelContinuum
end

function ProjectionTensor4(i::Int64, j::Int64, k::Int64, l::Int64)
    return (i == k ? 1.0 : 0.0)*(j == l ? 1.0 : 0.0)-1.0/3.0*(i == j ? 1.0 : 0.0)*(k == l ? 1.0 : 0.0)
end


function getProjectionTensor()
    mapDict = getTensorMapping()
    P = zeros(9,9)
    for k ∈ 1:3
      for l ∈ 1:3
          kl = getVoigtIndex(mapDict, k, l)
          for i ∈ 1:3
              for j ∈ 1:3
                  ij = getVoigtIndex(mapDict, i, j)
                  P[ij,kl] += ProjectionTensor4(i,j,k,l)
              end
          end
      end
  end
  return P
end

function get_σₘ_𝐬_mandel(σ_mandel::Array{T,1}) where T
    σₘ::Float64 = sum(σ_mandel[1:3])
    #Deviatoric Stress
    σ_mandel -= 1.0/3.0*σₘ*[1, 1, 1, 0, 0, 0, 0, 0, 0]
    𝐬::Float64 = sqrt(3.0/2.0)*norm(σ_mandel)
    return σₘ, 𝐬
end

function get_ϵₘ_𝒆_mandel(ϵ::Array{T,1}) where T
    ϵₘ::Float64 = sum(ϵ[1:3])
    #Deviatoric Stress
    ϵ -= 1.0/3.0*ϵₘ*[1, 1, 1, 0, 0, 0, 0, 0, 0]
    𝒆::Float64 = sqrt(2.0/3.0)*norm(ϵ)
    return ϵₘ, 𝒆
end

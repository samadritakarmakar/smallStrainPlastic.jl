# SmallStrainPlastic Documentation

SmallStrainPlastic is a library that aims to simplify the addition of different plastic models. The usage may be for a material point to judge the behaviour of a plastic model or to in a finite element model to simulate the behaviour of a plastic material.

The assumptions made while developing the library are the following:

The flow rule for plastic strain is considered to be:

𝛆̇ᵖ = λ̇  Θ(𝛔, 𝐪)

The flow rule for internal variable is considered to be:

𝛂̇ = λ̇  𝐡(𝛔, 𝐪)

The stress is considered as:

𝛔 = ℂ:(𝛆 - 𝛆ᵖ)

The hardening variable is considered as:

𝐪 = - 𝓗(𝛂)

Special Unicode characters using eg: "\bb*", "\bf*", or "\bsrc*" are used to define functions. It is highly recommended such unicode characters are avoided when defining internal variables to avoid confusion.

## Plasticity Model

```@docs
	SmallStrainPlastic.PlasticModel
```

## State of Plasticity

```@docs
	SmallStrainPlastic.State
```

```@docs
	SmallStrainPlastic.createStateDict
```

```@docs
	SmallStrainPlastic.updateStateDict!
```

```@docs
	SmallStrainPlastic.getState!
```
## Easier Finding of Jacobians

```@docs
	SmallStrainPlastic.denseJacobian!
```

```@docs
	SmallStrainPlastic.denseJacobian
```

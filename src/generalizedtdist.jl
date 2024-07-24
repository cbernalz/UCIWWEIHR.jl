# credit to Damon Bayer for all of these helper functions
"""
Create a generalized t distribution with location and scale.

# Arguments
- `μ`: Location parameter.
- `σ`: Scale parameter.
- `ν`: Degrees of freedom.

# Returns
A generalized t distribution object.

"""

function GeneralizedTDist(μ, σ, ν)
  if ν <= zero(ν)
    ν = nextfloat(zero(ν))
  end

  if σ <= zero(σ)
    σ = nextfloat(zero(σ))
  end

  Generalized_T = μ + σ*Distributions.TDist(ν)

  return(Generalized_T)

end

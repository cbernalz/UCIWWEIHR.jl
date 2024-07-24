# credit to Damon Bayer for all of these helper functions
"""
Create a re-parametrized negative binomial distribution in terms of mean and overdispersion.

# Arguments
- `μ`: Mean of the distribution.
- `ϕ`: Overdispersion parameter.

# Returns
A `Distributions.NegativeBinomial` distribution object.

"""
function NegativeBinomial2(μ, ϕ)
  p = 1 / (1 + μ / ϕ)

  if p <= zero(p)
    p = eps(zero(p))
  end

  if p >= one(p)
    p = prevfloat(one(p))
  end

  r = ϕ

  if r <= zero(r)
    r = eps(zero(r))
  end

  Distributions.NegativeBinomial(r, p)
end
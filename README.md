# NBConvolution

[![Build Status](https://travis-ci.org/slundberg/NBConvolution.jl.svg?branch=master)](https://travis-ci.org/slundberg/NBConvolution.jl)

This implements the distribution resulting from the convolution of several Negative binomial distributions (the distribution of the sum `S` in `S = X_1 + ... + X_n` when each `X_i` is a random variable following a Negative binomial distribution). The methods described in two papers are implemented. The first method by Furman represents an approximate solution. The second by Vallaisamy is an exact solution (that uses potentially large sums) and is the default option exported by this model:

"On the sums of compound negative binomial and gamma random variables" by Vallaisamy, 2009
http://projecteuclid.org/euclid.jap/1238592129

## Installation

```julia
Pkg.clone("https://github.com/slundberg/NBConvolution.jl.git")
```

## Usage

```julia
using NBConvolution

d1 = FurmanNegativeBinomialConvolution([4.0, 10.0], [0.5, 0.5], 1000)
d2 = NegativeBinomialConvolution([4, 10], [0.5, 0.5])
@assert abs(pmf(d1, 30) - pmf(d2, 30)) < 1e-10
```

You can also mix this with the `Distributions` package:

```julia
using Distributions

d2 = NegativeBinomialConvolution([4, 10], [0.5, 0.5])
dmix = MixtureModel([DiscreteUniform(0, 0), d2], [0.1, 0.9])
pdf(dmix, 10)
```


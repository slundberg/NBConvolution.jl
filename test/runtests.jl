using NBConvolution
using Distributions
using Base.Test

# compare the two methods to verify accuracy
d1 = FurmanNegativeBinomialConvolution([4.0, 1.0], [0.01, 0.5], 1000)
d2 = NegativeBinomialConvolution([4, 1], [0.01, 0.5])
@test abs(pdf(d1, 10) - pdf(d2, 10)) < 1e-10

d1 = FurmanNegativeBinomialConvolution([4.0, 10.0], [0.5, 0.5], 1000)
d2 = NegativeBinomialConvolution([4, 10], [0.5, 0.5])
@test abs(pdf(d1, 30) - pdf(d2, 30)) < 1e-10

d1 = FurmanNegativeBinomialConvolution([4.0, 10.0, 3.0, 6.0], [0.4, 0.5, 0.3, 0.6], 1000)
d2 = NegativeBinomialConvolution([4, 10, 3, 6], [0.4, 0.5, 0.3, 0.6])
@test abs(pdf(d1, 30) - pdf(d2, 30)) < 1e-10

# cdf
d2 = NegativeBinomialConvolution([4, 10], [0.7, 0.5])
@test abs(cdf(d2, 1000000) - 1.0) < 1e-10

# non-integer rs
d3 = NegativeBinomialConvolution([4.2, 10], [0.7, 0.5])
@test abs(cdf(d3, 1000000) - 1.0) < 1e-10
@test abs(pdf(d2, 10) - pdf(d3, 10)) > 1e-10
d4 = FurmanNegativeBinomialConvolution([4.2, 10.0], [0.7, 0.5], 1000)
@test abs(pdf(d3, 4) - pdf(d4, 4)) < 1e-10

# integrate with the distributions package
modelInner = MixtureModel([NegativeBinomialConvolution([3, 4], [0.1, 0.5]), NegativeBinomial(5, 0.4)], [0.6, 0.4])
pdf(modelInner, 0)
pdf(modelInner, 10)

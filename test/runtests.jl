using NegativeBinomialConvolution
using Base.Test

# compare the two methods to verify accuracy
d1 = FurmanNegativeBinomialSum([4.0, 1.0], [0.01, 0.5], 1000)
d2 = NegativeBinomialSum([4, 1], [0.01, 0.5])
@test abs(pmf(d1, 10) - pmf(d2, 10)) < 1e-10

d1 = FurmanNegativeBinomialSum([4.0, 10.0], [0.5, 0.5], 1000)
d2 = NegativeBinomialSum([4, 10], [0.5, 0.5])
@test abs(pmf(d1, 30) - pmf(d2, 30)) < 1e-10

d1 = FurmanNegativeBinomialSum([4.0, 10.0, 3.0, 6.0], [0.4, 0.5, 0.3, 0.6], 1000)
d2 = NegativeBinomialSum([4, 10, 3, 6], [0.4, 0.5, 0.3, 0.6])
@test abs(pmf(d1, 30) - pmf(d2, 30)) < 1e-10

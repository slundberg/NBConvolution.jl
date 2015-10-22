using Distributions
#using Distributions: check_args

export NegativeBinomialConvolution, pmf, cdf

immutable NegativeBinomialConvolution <: DiscreteUnivariateDistribution
    rs::Array{Int64,1}
    ps::Array{Float64,1}
    mean::Float64
    stddev::Float64

    """
    This is based on the paper:
    "On the sums of compound negative binomial and gamma random variables" by Vallaisamy, 2009
    http://projecteuclid.org/euclid.jap/1238592129
    """
    function NegativeBinomialConvolution(rs::Array{Int,1}, ps::Array{Float64,1})
        # @check_args(NegativeBinomialConvolution, minimum(rs) > zero(rs[1]))
        # @Distributions.check_args(NegativeBinomialConvolution, minimum(ps) > zero(ps[1]))
        # @Distributions.check_args(NegativeBinomialConvolution, Base.maximum(ps) < one(ps[1]))

        # just add the means and variances to get the convolution mean and variance
        meanVal = sum([(ps[i]*rs[i])/(1-ps[i]) for i in 1:length(ps)])
        varVal = sum([(ps[i]*rs[i])/(1-ps[i])^2 for i in 1:length(ps)])

        new(rs, ps, meanVal, sqrt(varVal))
    end
end

@Distributions.distr_support NegativeBinomialConvolution 0 Inf


params(d::NegativeBinomialConvolution) = (d.rs, d.ps)


"Call f on all possible ways n non-negative integers can sum to equal s."
function allsums_apply(f::Function, n::Int, s::Int)
    if n == 1
        f([s])
    elseif n == 2
        for i in 0:s
            f(Int64[i,s-i])
        end
    elseif n == 3
        for i in 0:s
            for j in 0:s
                if i + j > s break end
                f(Int64[i,j,s-i-j])
            end
        end
    elseif n == 4
        for i in 0:s
            for j in 0:s
                if i + j > s break end
                for k in 0:s
                    if i + j + k > s break end
                    f(Int64[i,j,k,s-i-j-k])
                end
            end
        end
    elseif n == 5
        for i in 0:s
            for j in 0:s
                if i + j > s break end
                for k in 0:s
                    if i + j + k > s break end
                    for l in 0:s
                        if i + j + k + l > s break end
                        f(Int64[i,j,k,l,s-i-j-k-l])
                    end
                end
            end
        end
    else
        error("n = $n is not yet implemented.")
    end
end

"Compute the exact PMF at the given point with a sum."
function Distributions.pdf(d::NegativeBinomialConvolution, s::Int)
    total = 0.0

    # short cut for very large values way past our distribution
    if d.mean + 20*d.stddev < s
        return 0.0
    end

    function inner_term(m)
        val = 1.0
        for j in 1:length(d.ps)
            try
                val *= binomial(d.rs[j] + m[j] - 1, m[j])*d.ps[j]^d.rs[j]*(1-d.ps[j])^m[j]
            catch
                val *= binomial(BigInt(d.rs[j] + m[j] - 1), BigInt(m[j]))*d.ps[j]^d.rs[j]*(1-d.ps[j])^m[j]
            end
        end
        total += val
    end
    allsums_apply(inner_term, length(d.ps), s)
    total
end

function Distributions.cdf(d::NegativeBinomialConvolution, s::Int)
    total = 0.0
    for i in 0:s
        total += pdf(d, i)
    end
    total
end

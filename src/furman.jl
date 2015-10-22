
export FurmanNegativeBinomialConvolution, pmf

"""
This is based on the paper:
"On the convolution of the negative binomial random variables" by Edward Furman, 2006
http://www.math.yorku.ca/~efurman/PDFs/NBconv.pdf

This is an approximation that is primarily superseded by an exact solution by Vallaisamy in 2009.
"""
type FurmanNegativeBinomialConvolution
    rs::Array{Float64,1}
    ps::Array{Float64,1}
    numTerms::Int64
    pmax::Float64
    R::Float64
    rsum::Int64
    xi::Array{Float64,1}
    delta::Array{Float64,1}
end
function FurmanNegativeBinomialConvolution(rs::Array{Float64,1}, ps::Array{Float64,1}, numTerms=500)
    pmax = maximum(ps)
    @assert pmax < 1

    R = 1.0
    for j in 1:length(ps)
        R *= (((1-ps[j])*pmax) / ((1-pmax)*ps[j]))^(-rs[j])
    end

    xi = zeros(numTerms)
    for i in 1:numTerms
        xi[i] = sum([(rs[j] * (1 - (1-pmax)*ps[j]/((1-ps[j])*pmax))^i)/i for j in 1:length(ps)])
    end

    delta = zeros(numTerms)
    delta[1] = 1
    for k in 1:numTerms-1
        delta[k+1] = sum([i*xi[i]*delta[k+1-i] for i in 1:k])/k
    end

    FurmanNegativeBinomialConvolution(float(rs), ps, numTerms, pmax, R, sum(rs), xi, delta)
end

"Approximate the PMF at the given point with a truncated sum."
function pmf(d::FurmanNegativeBinomialConvolution, s::Int)
    total = 0.0
    lastv = 0.0
    for k in 0:d.numTerms-1
        v = d.delta[k+1]*exp(lgamma(d.rsum + s + k) - lgamma(d.rsum + k) - lfact(s)) * d.pmax^(d.rsum+k)*(1-d.pmax)^s
        total += v

        # cheap check to make sure we don't have significant terms left to add
        if k == d.numTerms-1 && ((lastv < v && v > 1e-50) || v > 1e-10)
            println("$lastv $v")
            error("Sum expansion in NegativeBinomialConvolution left error! Consider using more terms.")
        end
        lastv = v
    end
    total * d.R
end

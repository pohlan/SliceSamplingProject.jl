using Distributions


function metropolis(θ0, logdensity; nsims=10^4, MH_width=0.75)
    θs   = zeros(nsims)
    acceptances = 0
    curθ = θ0
    for sim in 1:nsims
        propθ = rand(Uniform(curθ-MH_width, curθ+MH_width))
        accept = true
        if logdensity(propθ) < logdensity(curθ)
            R = exp(logdensity(propθ)-logdensity(curθ))
            coin_toss = rand(Binomial(1,R))
            if coin_toss == 0 accept = false end
        end
        if accept
            curθ = propθ
            acceptances += 1
        end
        θs[sim] = curθ
    end
    acceptance_ratio = acceptances / nsims * 100.
    return θs, acceptance_ratio
end


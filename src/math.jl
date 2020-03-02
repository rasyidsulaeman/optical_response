module math_module
using NumericalIntegration

export bisection, fermi, kramers_kronig

function bisection(f,a,b,itermax)
    for iter in 1:itermax
        if f(a)*f(b) > 0.
            println("Doesn't have a roots")
            break
        end

        c = (a+b)/2.

        if f(a)*f(c) < 0.
            b = c
        else
            a = c
        end

        if abs(b-a) < 0.00001
            return root = c
            println("convergent")
        end
    end
end

function fermi(w,μ,T)
    kT = 8.617332478e-5 * T
    if T == 0
        return float(w .< μ)
    else
        x = (w - μ) / kT
        return 1.0 / (1.0 + exp(x))
    end
end


function kramers_kronig(w,hw,fw)
    save_out = zeros(Float64,length(w))
    for (iw,ww) in enumerate(w)
        nums = ww^2 .- hw.^2
        select = filter(x -> !isapprox(x, 0, atol=1e-8),nums)
        if (length(select) == length(hw)-1)
            vals = (fw[2:end] .* ww) ./ select
            wintg = hw[2:end]
            save_out[iw] = integrate(wintg,vals) / π
        elseif (length(select) == length(hw))
            vals = (fw .* ww) ./ select
            wintg = hw
            save_out[iw] = 2*integrate(wintg,vals) / π
        end
    end
    return save_out
end

end

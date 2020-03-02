module TBmodel

# Importing the appropriate module
include("math.jl")
using .math_module

import Base.show

export k_generate, hamiltonian, tbmodel

#Generating kpoint mesh
function k_generate(nk)
    a = 1.0
    kk = collect(range(-π/a,length=nk,stop=π/a))
    kmesh = []
    for i in kk
        push!(kmesh,[i])
    end
    return kmesh
end

#Hamiltonian, it can be modify which fit with the problem in question
function hamiltonian(t,kmesh)
    H = zeros(Float64, length(kmesh))
    for (ik,k) in enumerate(kmesh)
        H[ik] = -2*t*(cos.(k[1]))
    end
    return H
end

#Tight Binding Model
function tbmodel(wrange::Tuple{Float64,Float64}, H, nfilling
                 ;nw = 501, T = 0.0, zeroplus = 0.012)

    @info "Calculating density of state"

    w = LinRange(wrange...,nw)

    dw = w[2] - w[1]

    function fmu(μ)
        intg = 0.0
        for iw in 1:length(w)
            if iw == 1 || iw == length(w)
                intg += fermi(w[iw],μ,T) * A0w[iw] * dw / 0.5
            else
                intg += fermi(w[iw],μ,T) * A0w[iw] * dw
            end
        end
        return nfilling - intg
    end

    Gw = zeros(ComplexF64, length(w))
    Gk = zeros(eltype(Gw), length(w), size(H)[1])
    for iw in 1:length(w)
        sumk = 0.0 + im*0.0
        for ik in 1:size(H)[1]
            dummy = 1.0 ./ ((w[iw] + im*zeroplus) - H[ik])
            Gk[iw,ik] = dummy
            sumk += dummy
        end
        Gw[iw] = sumk ./ size(H)[1]
    end

    A0w = -imag(Gw) ./ π

    μ = bisection(fmu,minimum(w),maximum(w),100)

    println("Chemical potential = ", μ)

    model = Dict("Green func w" => Gw,
                 "Green func k" => Gk,
                 "Chemical potential" => μ,
                 "Electron frequency" => w)
    return model
end

end

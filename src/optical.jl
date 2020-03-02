module Optical_Module

# Importing the appropriate module
include("math.jl")
using .math_module: kramers_kronig
using NumericalIntegration
using Interpolations

import Base.show

export vmatrices, dc_conductivity, complex_dielectric, loss_function, reflectance

function vmatrices(t,kmesh)
    vmat = zeros(Float64,length(kmesh))
    for (ik,kx) in enumerate(kmesh)
        vmat[ik] = 2*t*(sin.(kx[1]))
    end
    return vmat
end

function dc_conductivity(model::Dict,hwrange::Tuple{Float64,Float64},
                              vmat, kmesh; nhw = 201, T = 0.0)

    @info "Calculating Optical Conductivity"

    #Fundamental parameters
    e = 1.60217662e-19
    hbar = 1.0545718e-34
    a = 1.0e-9

    opt_prefact = π*e^2 / (hbar*a)

    hw = LinRange(hwrange...,nhw)

    #Importing some appropriate variable
    Gwk = model["Green func k"]
    μ = model["Chemical potential"]
    w = model["Electron frequency"]

    nw = length(w)
    dw = w[2] - w[1]

    rsigmaw = zeros(Float64,nhw)
    for (iω,ω) in enumerate(hw)

        #Integration boundary
        rangenu = (μ - ω,μ)
        meshnu = LinRange(rangenu...,nw)

        intg = zeros(Float64,length(meshnu))
        for (inu,nu) in enumerate(meshnu)

            #Index looping ν and ω in spectral function A0(k,ν) & A0(k,ν + ω)
            inu1 = floor(Int64,(nu - w[1])/dw) + 1
            inu2 = floor(Int64,(nu + ω - w[1])/dw) + 1

            if inu1 < 1 || inu2 > nw continue end
            sumk = 0.0
            for ik in 1:length(kmesh)
                Awk = -imag(Gwk[inu1,ik]) ./ π
                vA1 = vmat[ik] * Awk
                Awk = -imag(Gwk[inu2,ik]) ./ π
                vA2 = vmat[ik] * Awk
                vAvA = vA1 * vA2
                sumk += vAvA
            end

            if ω == 0.0
                dummy = sumk
            else
                dummy = sumk / ω
            end

            intg[inu] = dummy
        end

        vals = integrate(meshnu,intg, TrapezoidalEvenFast()) / length(kmesh)
        rsigmaw[iω] = vals
    end

    #Using interpolations to get the best part of integration
    #within Kramers-Kronig relations
    rsigma_itp = interpolate(rsigmaw,BSpline(Cubic(Line(OnGrid()))))
    sc_rsigma  = scale(rsigma_itp, hw)

    #Generate a new data point, which in this part will be 20*nhw data point
    hwnew      = LinRange(hwrange...,1002)
    rsigma_new = sc_rsigma(hwnew)
    isigmaw    = kramers_kronig(hw,hwnew,rsigma_new)

    #Get the optical conductivity
    sigmaw = opt_prefact .* (rsigmaw + im*isigmaw)

    @info "Calculation is done"

    opt_cond = Dict("photon_freq" => hw,
                    "conductivity" => sigmaw)

    return opt_cond
end

function complex_dielectric(opt_cond; epsinfty = 1.0)
    eps0 = 8.854187817e-12

    hw = opt_cond["photon_freq"]
    sigmaw = opt_cond["conductivity"]

    epsw = epsinfty .+ im*(sigmaw ./ (eps0 .* hw .* 1.5193e15))

    return epsw
end

function loss_function(epsw)
    epsw1 = real(epsw)
    epsw2 = imag(epsw)

    lw = epsw2 ./ (epsw1.^2 .+ epsw2.^2)
    return lw
end

function reflectance(epsw)
    kw = sqrt.(epsw)
    rw = abs.((kw .- 1).^2 ./ (kw .+ 1).^2)
    return rw
end

end

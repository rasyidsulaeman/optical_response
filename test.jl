"""
This is just an example how to implement tight-binding model along with
optical conductivity calculations
However, this code just work for T = 0 K
"""

#Importing the module
include("src/tb.jl")
include("src/optical.jl")

using .TBmodel
using .Optical_Module

using BenchmarkTools
using Plots; pyplot()


"""
Density of states calculations
Input : t => hopping parameter (eV)
        nfilling => number of electron in the band
        kmesh => mesh k-point
        H => hamiltonian
Output : Gw => Green function which it can be continue to calculate density of states as
         D0w = -imag(Gw) ./ π
"""

#Input some important parameters
t = 1.0
nfilling = 0.8

#First: we must generate k-point (which a point in brillouin zone)
#       and the hamiltonian
kmesh = k_generate(1000)
H     = hamiltonian(t,kmesh)

wrange = (-8.0,8.0)

#Create a model => output as a Dictionary
model = tbmodel(wrange,H,nfilling)

#Let's check our density of states
Gw = model["Green func w"]
w  = model["Electron frequency"]

D0w = -imag(Gw) ./ π

#Density of states plot
plot(w, D0w)
xlims!((wrange[1],wrange[2]))
ylims!((0.0,maximum(D0w)))
xlabel!("ω (eV)")
ylabel!("D (ω)")

"""
Let's try to calculate optical conductivity
Input : hwrange => range of photon frequency, typically between (0.0,4.0)
                   input must be a Tuple{T,T} where T <: Float64
        nhw => a number of photon frequency data, defaut : nhw = 201
        T => temperature, defaut : T = 0 K, unfortunately this code is
             intended for T = 0 K. Still in progress for T > 0 K
        model => A dictionary which include Green function with k-po*int grid
                 Gwk = G(k,ω)
Output : A dictionary which consists of
         hw => mesh of photon frequency
         σw => dc conductivity
"""

#Generate velocity matrices
vmat = vmatrices(t,kmesh)

#Range of photon frequency, continued to calculate conductivity
hwrange  = (0.0,4.0)
opt_cond = dc_conductivity(model,hwrange,vmat,kmesh)

"""
There is a parameters which imporant enough in calculation complex dielectric function,
which is ϵ1(∞), it can be any number as long as a float,
defaut value is 1.0, or ϵ1(∞) = 1.0

Ex: epsw = complex_dielectric(opt_cond, epsinfty = 1.0)
"""

#complex dielectric function
epsw = complex_dielectric(opt_cond)

#loss function
lw = loss_function(epsw)

#reflectance
rw = reflectance(epsw)

"""
Let's plots these variable
"""

hw = opt_cond["photon_freq"]
sigmaw = opt_cond["conductivity"]

#Conductivity
plot(hw, real(sigmaw) .* 1e-5, ylabel="Reσ (ω)") # we plot in unit of 10^3 Ohm cm^-1
plot(hw, imag(sigmaw) .* 1e-5, ylabel="Imσ (ω)")

#Complex dielectric function
plot(hw, real(epsw), ylabel="Reϵ (ω)")
plot(hw, imag(epsw), ylabel="Imϵ (ω)")

#Loss function
plot(hw, lw, ylabel="L (ω)")

#Reflectance
plot(hw, rw, ylabel="R (ω)")

xlims!((hwrange[1],hwrange[2]))
xlabel!("ω (eV)")

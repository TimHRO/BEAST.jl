using BEAST, CompScienceMeshes, SphericalScattering, LinearAlgebra, StaticArrays

f = 1e8
c = 2.99792458e8
μ = 4π * 1e-7
κ = 2π * f/c
spRadius = 1.0
r = 10.0
ϑ = range(0.0, stop=0.9999999999999999*π, length=18)  # 10° steps
ϕ = range(0.0, stop=0.9999999999999999*2π, length=36) # 10° steps

P = [SVector(cos(φ) * sin(θ), sin(φ) * sin(θ), cos(θ)) for θ in ϑ, φ in ϕ]
points_cartNF = P .* r
points_cartFF = P

# Solve the scatering problem using BEM

Γ = meshsphere(spRadius,0.2)
X = raviartthomas(Γ)

t = Maxwell3D.singlelayer(;wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
𝑒 = (n × E) × n

SL = Maxwell3D.singlelayer(; wavenumber=κ)

e = -assemble(𝑒, X)
T = assemble(SL, X, X)
u = T \ e

EF_MoM = potential(MWSingleLayerField3D(t), points_cartNF, u, X)
HF_MoM = potential(BEAST.MWDoubleLayerField3D(wavenumber=κ), points_cartNF, u, X) / (c * μ)
FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(t), points_cartFF, u, X)

# Solve the scattering problem by computing the Mie series

sp = PECSphere(radius=spRadius)
ex = planeWave(frequency=f)

EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
HF = scatteredfield(sp, ex, MagneticField(points_cartNF))
FF = scatteredfield(sp, ex, FarField(points_cartFF))

# Relative worst case errors

diff_EF = round(maximum(norm.(EF - EF_MoM) ./ maximum(norm.(EF))) * 100, digits=4)
diff_HF = round(maximum(norm.(HF - HF_MoM) ./ maximum(norm.(HF))) * 100, digits=4)
diff_FF = round(maximum(norm.(FF - FF_MoM) ./ maximum(norm.(FF))) * 100, digits=4)

print("E-field error: $diff_EF %\n")
print("H-field error: $diff_HF %\n")
print("FF-field error: $diff_FF %\n")
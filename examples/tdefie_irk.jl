using CompScienceMeshes, BEAST, StaticArrays, LinearAlgebra

radius = 1.0
Γ = meshsphere(radius=1.0, h=0.35)

∂Γ = boundary(Γ)

setminus(A,B) = submesh(!in(B), A)

edges = setminus(skeleton(Γ,1), ∂Γ)
verts = setminus(skeleton(Γ,0), skeleton(∂Γ,0))

Σ = Matrix(connectivity(Γ, edges, sign))
Λ = Matrix(connectivity(verts, edges, sign))

I = LinearAlgebra.I
PΣ = Σ * pinv(Σ'*Σ) * Σ'
PΛH = I - PΣ

X = raviartthomas(Γ)
sol = 1.0
Δt, Nt = 10.0, 200

(A, b, c) = butcher_tableau_radau_3stages()
T = StagedTimeStep(Δt, Nt, c, A, b, 5, 1.001)
V = X ⊗ T

Ip = diagm(@SVector ones(size(b,1)))

ℙΣ = kron(PΣ,Ip)
ℙΛH = kron(I - PΣ,Ip)

duration = 2 * 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
#Plots.plot(gaussian.(range(0,Nt*Δt,length=Nt)))

direction, polarisation = ẑ , x̂
E = planewave(polarisation, direction, BEAST.derive2(gaussian), sol)

T = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)

@hilbertspace j
@hilbertspace j′
tdefie_irk = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V
xefie_irk = solve(tdefie_irk)

Idp = diagm(ones(size(PΣ,1)))
Ã = inv(kron(Idp,A))
b̃ = kron(Idp,ones(size(b,1))*b')

j_sol = zeros(size(xefie_irk))
j_nonsol = zeros(size(xefie_irk))

a = 1.0

for i = 2:Nt
    j_sol[:,i] = ℙΛH * xefie_irk[:,i]
    j_nonsol[:,i] = a/(sol * Δt) * ℙΣ *(Ã*xefie_irk[:,i]-Ã*b̃*Ã*xefie_irk[:,i-1])
end

jf = j_sol[1:3:end,:] .+ j_nonsol[1:3:end,:]

import Plots, Plotly
Plots.plot((jf[1,:]))


Xefie_irk, Δω, ω0 = fouriertransform(jf, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[102]
print(ω1)
ue_irk = Xefie_irk[:,102] / fouriertransform(gaussian, numdiff=1)(ω1)

SL = Maxwell3D.singlelayer(; wavenumber=ω1/sol)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=ω1)
𝑒 = (n × E) × n
e = -assemble(𝑒, X)
T = assemble(SL, X, X)
ue = T \ e

fcr, geo = facecurrents(ue_irk, X)
Plotly.plot(patch(geo, norm.(fcr)))

Plots.plot(norm.(ue))

using SphericalScattering

c = 2.99792458e8

Θ = range(0, stop=2π, length=100)
Φ = 0.0
P = [ [cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)] for θ in Θ for ϕ in Φ]
points_cartNF = 5.0 .* P
points_cartFF = P

f = ω1/2π
μ = 4π * 1e-7

𝑇 = Maxwell3D.singlelayer(; wavenumber=ω1/sol)
EN_irk = potential(MWSingleLayerField3D(𝑇), points_cartNF, ue_irk, X)
EF_irk = -im * f/(2 * c) * potential(MWFarField3D(𝑇), points_cartFF, ue_irk, X)
EN_mom = potential(MWSingleLayerField3D(𝑇), points_cartNF, ue, X)
HF_mom = potential(BEAST.MWDoubleLayerField3D(wavenumber=ω1/sol), points_cartNF, ue, X) / (c * μ)
EF_mom = -im * f/(2 * c) * potential(MWFarField3D(𝑇), points_cartFF, ue, X)

sphere = PECSphere(radius=radius)
exc = planeWave(frequency=f)
EN_mie = scatteredfield(sphere, exc, ElectricField(points_cartNF))
HF_mie = scatteredfield(sphere, exc, MagneticField(points_cartNF))
EF_mie = scatteredfield(sphere, exc, FarField(points_cartFF))

diff_EF1 = maximum(norm.(EN_mie - EN_mom) ./ maximum(norm.(EN_mie)))
diff_EF2 = maximum(norm.(c.*EF_mie - EF_mom) ./ maximum(norm.(c.*EF_mie))) 
diff_HF = maximum(norm.(HF_mie - HF_mom) ./ maximum(norm.(HF_mie)))

fcr, geo = facecurrents(ue, X)
Plotly.plot(patch(geo, norm.(fcr)))

Plots.plot(Θ, real.(getindex.(EF_mom,1)), label="MoM")
Plots.plot!(Θ, real.(getindex.(c.*EF_mie,1)), label="MIE")

Plots.plot(norm.(jf[1,10:170]), yscale=:log10)
using CompScienceMeshes, BEAST, StaticArrays, LinearAlgebra

#compute mesh and projectors
radius = 1.0
Γ = meshsphere(radius=1.0, h=0.25)

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

#set speed of light to 1.0
sol = 1.0

#temporal step size and original number of time steps (we are actually solving for p*Nt time steps)
Δt, Nt = 1.0, 400

#choose suitable Runge-Kutta tableau
(A, b, c) = butcher_tableau_radau_3stages()

#struct that contains all temporal information
T = StagedTimeStep(Δt, Nt, c, A, b, 10, 1.0001)

#space-time basis that combines temporal and spacial information
V = X ⊗ T

Ip = diagm(@SVector ones(size(b,1)))

#blow up projectors according to the number of stages
ℙΣ = kron(PΣ,Ip)
ℙΛH = kron(I - PΣ,Ip)

#create gaussian pulse
duration = 2 * 40 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
Plots.plot(gaussian.(range(0,Nt*Δt,length=Nt)))

#gaussian pulse defines aplitude of planewave
direction, polarisation = ẑ , x̂ 
E = planewave(polarisation, direction, BEAST.derive(gaussian), sol)

#use differentiated time-domain single layer operator
T = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)

@hilbertspace j
@hilbertspace j′
tdefie_irk = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V

#solve for the auxiliary current yᵢ
xefie_irk = solve(tdefie_irk)

#retreive originial current jᵢ
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

#get rid of intermediate stages
jf = j_sol[1:3:end,:] .+ j_nonsol[1:3:end,:]

import Plots, Plotly
Plots.plot((j_sol[1,:]))
Plots.plot!((j_nonsol[1,:]))

Xefie_irk, Δω, ω0 = fouriertransform(jf, Δt, 0.0, 2)
J_nonsol, _, _ = fouriertransform(j_nonsol[1:3:end,:], Δt, 0.0, 2)
J_sol, _, _ = fouriertransform(j_sol[1:3:end,:], Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmax(abs.(Xefie_irk[1,:]))

ω1 = ω[i1+1]
print(ω1)
ue_jf = Xefie_irk[:,i1+1] / fouriertransform(gaussian, numdiff=0)(ω1)
ue_nonsol = J_nonsol[:,i1+1] / fouriertransform(gaussian, numdiff=0)(ω1)
ue_sol = J_sol[:,i1+1] / fouriertransform(gaussian, numdiff=0)(ω1)

#Space Time Galerkin solution
T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)

tdefie = @discretise SL[j′,j] == -1.0E[j′]   j∈V  j′∈W
xefie = BEAST.motsolve(tdefie)

#Mie Series 
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
EN_irk = potential(MWSingleLayerField3D(𝑇), points_cartNF, ue_nonsol, X)
HF_irk = potential(BEAST.MWDoubleLayerField3D(wavenumber=ω1/sol), points_cartNF, ue_sol, X) / (c * μ)
EF_irk_nonsol = potential(MWFarField3D(𝑇), points_cartFF, ue_nonsol, X)
EF_irk_sol = potential(MWFarField3D(𝑇), points_cartFF, ue_sol, X)
EF_irk = -im * f/(2 * c) * (EF_irk_nonsol .+ EF_irk_sol)

sphere = PECSphere(radius=radius)
exc = planeWave(frequency=f)

EN_mie = scatteredfield(sphere, exc, ElectricField(points_cartNF))
HF_mie = scatteredfield(sphere, exc, MagneticField(points_cartNF))
EF_mie = scatteredfield(sphere, exc, FarField(points_cartFF))

diff_EN = maximum(norm.(EN_mie - EN_irk) ./ maximum(norm.(EN_mie)))
diff_EF = maximum(norm.(c.*EF_mie - EF_irk) ./ maximum(norm.(c.*EF_mie))) 
diff_HF = maximum(norm.(HF_mie - HF_irk) ./ maximum(norm.(HF_mie)))

diff_MoM = maximum(norm.(xefie[:,1:200]-jf[:,1:200]) ./ maximum(norm.(xefie[:,1:200])))

fcr1, geo1 = facecurrents(j_nonsol[1:3:end,120], X)
fcr2, geo2 = facecurrents(j_sol[1:3:end,120], X)

p1 = Plotly.plot(patch(geo1, norm.(fcr1)))
p2 = Plotly.plot(patch(geo2, norm.(fcr2)))
p = [p1 p2]

ti = range(0,Nt*Δt,length=Nt)

Plots.plot(Θ, real.(getindex.(EF_irk,1)), label="IRK")
Plots.plot!(Θ, real.(getindex.(c.*EF_mie,1)), label="MIE")

Plots.plot(ti = range(0,Nt*Δt,length=Nt), norm.(jf[1,:]), xlabel = "Time t in s", ylabel = "Surface current densitiy |j(t)| in A/m", label = "IRK")
Plots.plot!(norm.(xefie[1,:]), label = "MoT")

Plots.plot(ti[3:390], norm.(jf[1,3:390]), xlabel = "Time t in s", ylabel = "Surface current densitiy |j(t)| in A/m", yscale=:log10, label = "IRK")
Plots.plot!(norm.(xefie[1,3:390]), yscale=:log10, label = "MoT")





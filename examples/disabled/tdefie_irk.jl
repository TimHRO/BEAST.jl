using CompScienceMeshes, BEAST, StaticArrays, LinearAlgebra, Plots

Γ = meshsphere(radius=1.0, h=0.35)
X = raviartthomas(Γ)
sol = 1.0
Δt, Nt = 10.0, 200

(A, b, c) = butcher_tableau_radau_2stages()
T = StagedTimeStep(Δt, Nt, c, A, b, 10, 1.001)
V = X ⊗ T

duration, delay, amplitude = 2 * 20 * Δt, 2 * 30 * Δt, 1.0
gaussian = creategaussian(duration, delay, amplitude)

direction, polarisation = ẑ , x̂
E = planewave(polarisation, direction, BEAST.derive(gaussian), sol)

T = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)

@hilbertspace j
@hilbertspace j′
tdefie_irk = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈V
xefie_irk = solve(tdefie_irk)

j = xefie_irk[1:2:end,:]

import Plots
Plots.plot(j[1,:])

import Plotly
fcr, geo = facecurrents(j[:,10], X)
Plotly.plot(patch(geo, norm.(fcr)))


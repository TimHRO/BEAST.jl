using CompScienceMeshes, BEAST, LinearAlgebra, Plots

Γ = meshsphere(radius=1.0, h=0.35)
X = raviartthomas(Γ)

sol = 1.0
Δt = 10000.0
Nt = 200

duration = 2 * 70 * Δt
delay =  0.7 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
y = derive(gaussian).(range(0,Nt*Δt,length=Nt))
Plots.plot(range(0,Nt*Δt,length=Nt),y)
direction, polarisation = ẑ, x̂

E = planewave(polarisation, direction, BEAST.derive2(gaussian), sol)

#tdefie---------------------------------------

T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

@hilbertspace j
@hilbertspace j′

SL = TDMaxwell3D.singlelayer(speedoflight=sol, numdiffs=1)

tdefie = @discretise SL[j′,j] == -1.0E[j′]   j∈V  j′∈W
xefie = solve(tdefie)
Plots.plot(xefie[2,:])

Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
#ω = collect(ω0 .+ (0:Nt-1)*Δω)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
ωp = ω.+ω[end]
ffdg = fouriertransform(gaussian, numdiff=1).(ω)
Plots.plot(ω, real(ffdg))


#efie solution for ω ----------------------------

f_efie = zeros(ComplexF64, size(xefie))

for (i,w) in enumerate(ω)
    if w == 0.0
        f_efie[:,i] .= 0.0
        continue
    end
    κ = w/sol
    print(κ)
    t = Maxwell3D.singlelayer(;wavenumber=κ)
    E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
    𝑒 = (n × E) × n

    SL = Maxwell3D.singlelayer(; wavenumber=κ)

    e = -assemble(𝑒, X)
    T = assemble(SL, X, X)
    f_efie[:,i] = (T \ e) * fouriertransform(gaussian, numdiff=1)(κ)
end

#inverse fouriertransform and demodulate

inverse_efie, Δti, t0i = BEAST.inversefouriertransform(f_efie, Δω, 0.0, 2)
t = collect(0.0.+ (0:Nt-1)*Δti)
e = exp.(-im*ω[begin]*t)
Plots.plot!(real(inverse_efie[1,:]))
fefie,_,_ = fouriertransform(xefie_irk, Δt, 0.0, 2)
Plots.plot(xefie_irk[1,:])
td_efie = zeros(size(inverse_efie))
for i in 1:size(xefie,1)
    td_efie[i,:] = real(inverse_efie[i,:] .*e)
end

import Plots, Plotly

diff_MOT_max = norm((norm.(xefie - td_efie))) ./ norm(norm.(td_efie))
diff_IRK_max = norm((norm.(jf - td_efie))) ./ norm(norm.(jf))

Plots.plot(t,td_efie[1,:], label = "efie solution")
Plots.plot!(t,xefie[1,:], label = "MOT solution")
Plots.plot!(t,jf[1,:], label = "IRK solution")
xlabel!("t")
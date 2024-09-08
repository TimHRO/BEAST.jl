using Plots
using JLD2

eigenval = load("eigenvalues.jld2", "ev")
function my_circle(r,h,w)
    ϕ = range(0,2π, 1000)
    x = w .+ r.*cos.(ϕ)
    y = h .+ r.*sin.(ϕ)
    return x,y
end

Plots.plot(my_circle(1.0,0,0))
Plots.scatter!(eigenval)
print(maximum(norm.(eigenval))/minimum(norm.(eigenval)))
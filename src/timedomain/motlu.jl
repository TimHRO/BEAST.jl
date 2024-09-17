#JLD2 is used for saving eigenvalues of Zc - computation of eigenvalues is not neccessary for solution 
#using JLD2
motsolve(eq) = td_solve(eq)

function td_solve(eq)

    V = eq.trial_space_dict[1]

    A = assemble(eq.equation.lhs, eq.test_space_dict, eq.trial_space_dict)
    T = eltype(A)
    S = zeros(T, size(A)[1:2])
    ConvolutionOperators.timeslice!(S, A, 1)

    iS = inv(S)
    b = assemble(eq.equation.rhs, eq.test_space_dict)

    #Compute eigenvalues of reccurence relation
    #Dély, A., F.P. Andriulli, and K. Cools. 2018. “Stable TD-EFIE Discretized with Implicit Runge-Kutta Methods.”

    time_info = temporalbasis(V)
    #Nconv =  time_info.zTransformedTermCount
    #sA = size(A,1)
    #Ac = zeros(sA*Nconv, sA*Nconv)
    #Zc = zeros(sA*Nconv, sA*Nconv)
    #for i in 1:Nconv-1
    #    Zc[1:sA, 1+(i-1)*sA:i*sA]=-iS*ConvolutionOperators.timeslice(A,i+1)
    #end
    #for i in 1:Nconv-2
    #    Zc[1+sA+(i-1)*sA:(i+1)*sA, 1+(i-1)*sA:i*sA]=Matrix(I,sA,sA)
    #end
    #ev = eigen(Zc).values
    #@save "eigenvalues.jld2" ev
    nt = numfunctions(temporalbasis(V))
    marchonintime(iS, A, b, nt)
end

"""
    marchonintime(W0,Z,B,I; convhist=false)

Solve by marching-on-in-time the causal convolution problem defined by `(W0,Z,B)`
up to timestep `I`. Here, `Z` is an array of order 3 that contains a discretisation
of a time translation invariant retarded potential operator. `W0` is the inverse of
the slice `Z[:,:,1]`.

Keyword arguments:
    - 'convhist': when true, return in addition to the space-time data for the
    solution also the vector of convergence histories as returned each time step
    by the supplied solver `W0`.
"""
function marchonintime(W0,Z,B,I; convhist=false)

    T = eltype(W0)
    M,N = size(W0)
    @assert M == size(B,1)

    x = zeros(T,N,I)
    y = zeros(T,N)
    csx = zeros(T,N,I)

    ch = []
    for i in 1:I
        R = B[:,i]
        k_start = 2
        k_stop = I

        fill!(y,0)
        ConvolutionOperators.convolve!(y,Z,x,csx,i,k_start,k_stop)
        b = R - y
        xi, chi = BEAST.solve(W0, b)
        x[:,i] .+= xi
        push!(ch, chi)
        if i > 1
            csx[:,i] .= csx[:,i-1] .+ x[:,i]
        else
            csx[:,i] .= x[:,i]
        end

        (i % 10 == 0) && print(i, "[", I, "] - ")
    end

    if convhist
        return x, ch
    else
        return x
    end
end

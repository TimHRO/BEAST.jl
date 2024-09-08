

struct Gaussian{T}
    scaling::T
    width::T
    delay::T
end

struct TayloredGaussian{T}
    scaling::T
    width::T
    delay::T
end

Gaussian(;scaling=1.0, width, delay) = Gaussian(typeof(width)(scaling), width, delay)
function (g::Gaussian)(s::Real) 
    y = 4*g.scaling/(g.width*√π) * exp(-(4*(s-g.delay)/g.width)^2) #* cos(2π*1e5*s)
    if exp(-(4*(s-g.delay)/g.width)^2) == 1.0
        print(s-g.delay)
        print(" ")
        #print("maschine precision")
    end
    return y
end

TayloredGaussian(;scaling=1.0, width, delay) = TayloredGaussian(typeof(width)(scaling), width, delay)
function (g::TayloredGaussian)(s::Real)
    return TaylorSum(g,s,30)
end

function TaylorSum(g::TayloredGaussian,s,n)
    y = -(4*(s-g.delay)/g.width)^2
    if s-g.delay == 0.0
        print("maschine precision in taylor")
    end
    mysum  = 0
    for i in 0:n
        mysum += y^i/factorial(big(i)) 
    end
    return 4*g.scaling/(g.width*sqrt(π)) * mysum
end


function creategaussian(width,s0,scaling=one(typeof(width)))
    #f(s) = 4*scaling/(width*sqrt(π)) * exp(-(4*(s-s0)/width)^2)
    Gaussian(scaling, width, s0)
    #f(s) = scaling * exp(-(4*(s-s0)/width)^2)
end

function createTayloredGaussian(width,s0,scaling=one(typeof(width)))
    TayloredGaussian(scaling, width, s0)
end


function  fouriertransform(g::Gaussian; numdiff=0)
    scaling = g.scaling
    width = g.width
    s0 = g.delay
    ft(w) = (im*w)^numdiff * scaling * exp(-im*w*s0 - (width*w/8)^2) / sqrt(2π)
end


function fouriertransform(a::Array, dt, t0, dim=1)
    n = size(a,dim)
    dω = 2π / (n*dt)
    b = fftshift(fft(a, dim), dim) * dt / sqrt(2π)
    ω0 = -dω * div(n,2)
    b, dω, ω0
end

function inversefouriertransform(a::Array, dω, ω0, dim=1)
    n = size(a,dim)
    dt = 2π/ (n*dω)
    b = ifft(a,dim) * sqrt(2π) / dt
    t0 = -dt * div(n,2)
    b, dt, t0
end

fouriertransform(a::Array; stepsize, offset, dim=1) = fouriertransform(a, stepsize, offset, dim)


derive(g::Gaussian) =  s -> g(s) * (-8 * (s-g.delay)/g.width) * (4/g.width)

derive2(g::Gaussian) = s -> -(8*4)/g.width^2 * (g(s) + (s-g.delay) * derive(g)(s))


struct ErrorFunction{T}
    scaling::T
    width::T
    delay::T
end

function (f::ErrorFunction)(s)
    #y = f.scaling * 0.5 * (1 + erf(4*(s-f.delay)/f.width))
    y = 4*(s-f.delay)/f.width
    #sum = 0
    #for (i,k) in enumerate(1:2:200)
        #sum += (-1)^(i+1)*y^k/(k*factorial(big(i)))
        #print(i)
        #print(k)
        #print("  ")
    #end
    #if y == 0.5
    #    print(4*(s-f.delay)/f.width)
    #    print("  ")
    #end
    return f.scaling * 0.5 * (1 + 1 - 1/((1+0.278393*y+0.230389*y^2+0.000972*y^3+0.078108*y^4)^4))
    #return f.scaling * 0.5 * (1 + 2/sqrt(π) * sum)
    #return y
end

function integrate(f::Gaussian)
    return ErrorFunction(f.scaling, f.width, f.delay)
end

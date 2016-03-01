function chi2(X)
    #=
    Calculate the Chi2 value
    X[1,:] => Obs
    X[2,:] => Err
    X[3,:] => Model
    =#
    return sum(((X[1,:] - X[3,:]) / X[2,:])^2.)
end

function Merit(X)
    #= Given a Chi2 value
    we calculate a likelihood as our merit function
    =#
    return exp(-chi2(X)/2.)
end

function FuncLinear(x,Par)
  # linear function y = ax + b
  return Par[1].*x + Par[2]
end

function GenFakeData(N,F,P,E,mean_noise)
  #=
  N => Number of points
  F => Function used
  P => Paramters
  E => Average error on y-axis points
  mean_noise => Mean value of noise (Typically 0)
  =#
  noise = randn(N).*E
  x     = linspace(0,N,N)
  y     = F(x,P)+noise
  yerr  = ones(length(x))#.*E
  return x,y,yerr,FuncLinear(x,P)
end

function MCMC(X,F,P,Const,S,C)
    #=
    X => Data (y,yerr,model)
    F => Function used
    P => Parameters
    S => Scale
    C => Chain length
    =#

    L = Merit(X)
    moves = 0
    L_chain = zeros(length(P),C)

    for i = 1:C
        println("hello")
        if i%10 == 0.
            println( (i/C)*100.)#," % done"
        end
        jump = randn(1,length(S)).*S    # element-wise product of vectors .*
        P = P + jump
        new_fit =    1.     
    end
return 1.
end
# Initial Parameters
# Linear Function Example
Pin     = [0.018 .0]
Pfake   = [0.015 0.01]
step    = [0.001 0.016]
F       = FuncLinear

X = GenFakeData(100,F,Pfake,0.1,0.)

Pkg.add("PyPlot")
using PyPlot
plot(X[1], X[2], color="red")
show()
#println(MCMC(X,F,P,Const,S,C))

function chi2(X)
    #=
    Calculate the Chi2 value
    X[1,:] => x-axis
    X[2,:] => Obs
    X[3,:] => Err
    X[4,:] => Model
    =#
    return sum(((X[2,:] - X[4,:]) / X[3,:])^2.)

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

function GenFakeData(N,F,P,E)
  #=
  N => Number of points
  F => Function used
  P => Paramters
  E => Average error on y-axis points
  =#
  noise = randn(N).*E
  x     = linspace(0,N-1,N)
  y     = F(x,P)+noise
  yerr  = ones(length(x)).*E
  return [x y yerr FuncLinear(x,P)]'
end

function Distrib(x)
   # Finds median and 68% interval of array x. 
   y    =   sort(x)
   up   =   y[int(0.8413*length(y))]
   down =   y[int(0.1587*length(y))]
   med  =   y[int(0.5*length(y))]
   return [med up down]   
end

function PrintUncertainties(chain,P,S)
  # Prints the median parameter values
  # with their associated uncertainties
  println("====================================")
  for i = 1:length(P)
    if S[i] != 0
      println("\t\t+$(round(Distrib(chain[:,P[i]])[1],5))")
      println(round(Distrib(chain[:,P[i]])[1],5))
      println("\t\t-$(round(Distrib(chain[:,P[i]])[3],5))")
    end
  end
  println("====================================")
end

function PlotChain(chain,P,S)

  x = chain[:,P[1]]
  y = chain[:,P[2]]
  
  fig = figure()

  hist2d(x=x,y=y)
  show()
  #=
  # Top plot
  #top = subplot2grid((3,3), (0, 0), colspan=2)
  hist(chain[:,P[1]],bins=30)
  #axvline(Distrib(chain[:,P[0]])[0],color="red",lw=2)
  #axvline(Distrib(chain[:,P[0]])[1],color="red",lw=2,linestyle='--')
  #axvline(Distrib(chain[:,P[0]])[2],color="red",lw=2,linestyle='--')
  #top.get_xaxis().set_ticklabels([])
  minorticks_on()


  # Right hand side plot
  right = plt.subplot2grid((3,3), (1, 2), rowspan=2)
  right.hist(chain[:,P[1]],orientation='horizontal',bins=30)
  plt.axhline(Distrib(chain[:,P[1]])[0],color="red",lw=2)
  plt.axhline(Distrib(chain[:,P[1]])[1],color="red",lw=2,linestyle='--')
  plt.axhline(Distrib(chain[:,P[1]])[2],color="red",lw=2,linestyle='--')
  right.get_yaxis().set_ticklabels([])
  right.xaxis.set_major_locator(LinearLocator(5))
  plt.minorticks_on()
  
  # Center plot
  center = plt.subplot2grid((3,3), (1, 0), rowspan=2, colspan=2)
  center.hist2d(chain[:,P[0]],chain[:,P[1]],bins=30)
  plt.minorticks_on()
  
  # Corner plot
  corner = plt.subplot2grid((3,3), (0, 2))
  corner.get_xaxis().set_ticklabels([])
  corner.get_yaxis().set_ticklabels([])
  corner.plot(chain[:,P[0]],chain[:,P[1]],'-k')
  plt.minorticks_on()
  =#
end

function MCMC(X,F,P,S,C)
    #=
    X => Data (x, y,yerr,model)
    F => Function used
    P => Parameters
    S => Scale
    C => Chain length
    =#

    L       = Merit(X)
    moves   = 0
    chain   = zeros(C,length(P))
    LChain  = zeros(C,1)
    
    for i = 1:C
        if i%(C/10) == 0.
            println("$((i/C)*100.) % done")#," % done"
        end
        jump        = randn(1,length(S)).*S    # element-wise product of vectors .*
        P           = P + jump
        X[4,:]      = FuncLinear(X[1,:],P)
        Lnew        = Merit(X)     
        LChain[i]   = Lnew
        ratio       = Lnew/L
        if rand()   > ratio
            P       = P - jump
            moved   = 0
        else
            L       = Lnew
            moved   = 1
        end
        moves += moved
        chain[i,:] = P
    end
    println("\nAccepted steps: $(100.*(moves/C)) %")
    return chain, moves
end



# Initial Parameters
# Linear Function Example
Pin     = [0.018 0.0]
Pfake   = [0.015 0.01]
step    = [0.01 0.016]
F       = FuncLinear
C       = 100   # Number of iterations type: int


X       = GenFakeData(100,F,Pfake,0.1)

chain, moves   = MCMC(X,F,Pin,step,C)
Pout    = chain[moves,:]
#println(Pout)

P_plot = [1 2]
PrintUncertainties(chain,P_plot,step)
#using Plotly
#PlotChain(chain,P_plot,step)



#x = X[1,:]
#y = X[2,:]



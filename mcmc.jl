#using PyPlot
import PyPlot
const plt = PyPlot

function chi2(D,model)
    #=
    Calculate the Chi2 value
    D[1] => x-values
    D[2] => Obs
    D[3] => Err
    =#
    return sum(  ((  (D[2]-model) - zeros(length(D[1]))) / D[3])^2.  )
end

function Merit(D,model)
    #= Given a Chi2 value
    we calculate a likelihood as our merit function
    =#
    return exp(-chi2(D,model)/2.)
end

function FuncLinear(x,Par)
  # linear function y = ax + b
  return Par[1].*x + Par[2]
end

function GenData(N,F,P,E)
  #=
  N => Number of points
  F => Function used
  P => Paramters
  E => Average error on y-axis points
  =#
  println("Generating data...")
  noise = randn(N).*E
  x     = linspace(0,N-1,N)
  y     = F(x,P)+noise
  yerr  = ones(length(x)).*E
  return x,y,yerr
end

function Distrib(x)
   # Finds median and 68% interval of array x. 
   y    =   sort(x)
   up   =   y[round(Int,0.8413*length(y))]
   down =   y[round(Int,0.1587*length(y))]
   return median(y), up, down   
end

function PrintUncertainties(chain,P,S)
  # Prints the median parameter values
  # with their associated uncertainties

  A = zeros(length(P),3)
  
  for i = 1:length(P)
    A[i,3] = Distrib(chain[:,P[i]])[2]-Distrib(chain[:,P[i]])[1]
    A[i,1] = Distrib(chain[:,P[i]])[1]
    A[i,2] = Distrib(chain[:,P[i]])[1]-Distrib(chain[:,P[i]])[3]    
  end
  
  return A
end

function PlotFunction(x,y,err,Finit,Ffinal)
    plt.plot(x,Ffinal,color="red")
    plt.errorbar(x,y,yerr=err,linestyle ="None",ecolor="black")#,fmt='o',ecolor="black")
    plt.plot(x,y, marker="o",linestyle ="None", color="black")
    plt.plot(x,Finit,color="blue")
    plt.show()
end

function PlotChain(chain,P,S)
  println("Plotting...")
  fig = plt.figure()
  plt.clf()
  

  # Top plot
  top = plt.subplot2grid((3,3), (0, 0), colspan=2)
  plt.plt[:hist](chain[:,P[1]],bins=30)
  plt.axvline(Distrib(chain[:,P[1]])[1],color="red",lw=2)
  plt.axvline(Distrib(chain[:,P[1]])[2],color="red",lw=2,linestyle="--")
  plt.axvline(Distrib(chain[:,P[1]])[3],color="red",lw=2,linestyle="--")
  #plt.get_xaxis().set_ticklabels([])
  #plt.[:xaxis][:set_ticklabels]() # https://gist.github.com/gizmaa/7214002
  plt.minorticks_on()

  
  # Right hand side plot
  right = plt.subplot2grid((3,3), (1, 2), rowspan=2)
  plt.plt[:hist](chain[:,P[2]],orientation="horizontal",bins=30)
  plt.axhline(Distrib(chain[:,P[2]])[1],color="red",lw=2)
  plt.axhline(Distrib(chain[:,P[2]])[2],color="red",lw=2,linestyle="--")
  plt.axhline(Distrib(chain[:,P[2]])[3],color="red",lw=2,linestyle="--")
  #right.get_yaxis().set_ticklabels([])
  #right.xaxis.set_major_locator(LinearLocator(5))
  plt.minorticks_on()

  # Center plot
  center = plt.subplot2grid((3,3), (1, 0), rowspan=2, colspan=2)
  plt.hist2D(chain[:,P[1]],chain[:,P[2]],bins=30)
  plt.minorticks_on()
  plt.xlabel("a")
  plt.ylabel("b")
  
  # Corner plot
  corner = plt.subplot2grid((3,3), (0, 2))
  #corner.get_xaxis().set_ticklabels([])
  #corner.get_yaxis().set_ticklabels([])
  plt.plot(chain[:,P[1]],chain[:,P[2]],"-k")
  plt.minorticks_on()

  plt.show()
  #plt.savefig(Target+'_param_'+str(P[0])+'.pdf',paper='a4',orientation='landscape',dpi=300,transparent=True, bbox_inches='tight', pad_inches=0.1)

#=
x = chain[:,P[1]]
y = chain[:,P[2]]

data = [
  Dict(
    "x" => x,
    "y" => y,
    "type" => "histogram2d"
  )
]
response = Plotly.plot(data, Dict("filename" => "2d-histogram", "fileopt" => "overwrite"))
plot_url = response["url"]
=#
end

function MCMC(D,F,P,S,C)
    #=
    D => Data (x = D[1], y, yerr)
    F => Function used
    P => Parameters
    S => Step size
    C => Chain length
    =#
    println("Performing MCMC...")
    model   = F(D[1],Pin)
    L       = Merit(D,model)
    moves   = 0
    chain   = zeros(C,length(P))
    LChain  = zeros(C,1)

    for i = 1:C
        if i%(C/10) == 0.
            println("$((i/C)*100.) % done")#," % done"
        end
        jump        = randn(1,length(S)).*S    # element-wise product of vectors .*
        P           = P + jump
        NewModel    = F(D[1],P)
        Lnew        = Merit(D,NewModel)     
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
Pin     = [0.1 -1.0]
Pfake   = [0.1 -1.0] # Actual answer
step    = [0.0025 0.0025]

F       = FuncLinear
C       = 100000   # Number of iterations type: int


D       = GenData(50,F,Pfake,0.3)

#PlotFunction(D[1],D[2],D[3],F(D[1],Pin),F(D[1],Pfake))


chain, moves   = MCMC(D,F,Pin,step,C)
#Pout    = chain[moves,:]

Pplot = [1 2]

A       = PrintUncertainties(chain,Pplot,step,)

Pmedian = [A[1,1] A[2,1]]

PlotFunction(D[1],D[2],D[3],F(D[1],Pmedian),F(D[1],Pfake))

println("====================================")
println("\na\t=\t", round(A[1,1],4),"\t+",round(A[1,3],4),"\t-",round(A[1,2],4))
println("b\t=\t",     round(A[2,1],4),"\t+",round(A[2,3],4),"\t-",round(A[2,2],4))
println("====================================\n")

PlotChain(chain,Pplot,step)


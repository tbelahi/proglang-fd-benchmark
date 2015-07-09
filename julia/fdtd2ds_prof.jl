function make_source(source,fc,nt,dt)
  coeff=4
  tsour=1/fc
  t = [0:nt-1]*dt
  lsour=tsour/dt
  t0=tsour*1.5
  #T0=tsour*1.5
  tau = zeros(nt)::Array{Float64,1}
   for l=1:nt
    tau[l]=pi*(t[l]-t0)/t0
  end
  #tau=pi*(t-t0)/T0
  fs = zeros(nt)::Array{Float64,1}
  if source==0
     for k=1:nt
      fs[k] =(1-coeff*tau[k])*exp(-2*tau[k]*tau[k])
    end
  else
    fs=source
  end
  return fs
end

function fdtd2ds(source)
  # fdtd2ds: finitite differences time domain 2ds

  # Define velocity model
  const nx = 400
  const ny = 400
  const nx2 = nx+1
  const ny2 = ny+1
  const vp0 = 2000
  vp = vp0*ones(nx2,ny2)

  #vp = vp .* vp
   for j=1:ny2
     for i=1:nx2
       vp[i,j]=vp[i,j]*vp[i,j]
    end
  end
  const rho = 1500
  # Define source position and receiver position
  const isx = 200
  const isy = 100
  const irx = 200
  const iry = 300

  # source and fd parameters fd= finite difference
  const nt = 3000 #better is set at 3000
  const fc = 1000
  dx = vp0/(30*fc)
  dy = dx
  dt = dx/vp0*1/(sqrt(2)*2*pi)

  # Define source-time function
  fs = make_source(source,fc,nt,dt)

  # Define PML (perfectly matched layer boundaries) properties
  npml = 30
  pmlfac = 50
  pmlexp = 2
  qx = zeros(nx2,ny2)::Array{Float64,2}# inutile
  qy = zeros(nx2,ny2)::Array{Float64,2}# inutile
  qx = zeros(nx,ny)::Array{Float64,2}
  qy = zeros(nx,ny)::Array{Float64,2}
   for i=1:ny
     for a = 1:npml
       qx[a,i] = pmlfac*(npml-a)^pmlexp # left explained by the transposition of the matrices when plotted  x is vertical axis in the code and y is horizontal, the trickt to change it is to flip the matrices when plotting
       qx[nx-a+1,i] = pmlfac*(npml-a)^pmlexp # right
    end
  end
   for a=1:npml
     for i=1:nx
       qy[i,a] = pmlfac*(npml-a)^pmlexp # thp
       qy[i,ny-a+1] = pmlfac*(npml-a)^pmlexp # bottom
    end
  end

  qx = [qx[:,1] qx] # pour que ça fasse la taille nx2, ny2
  qx = [qx[1,:]; qx]
  qy = [qy[:,1] qy]
  qy = [qy[1,:]; qy]



  # Initialize fields
  px = zeros(nx2,ny2)::Array{Float64,2}
  py = zeros(nx2,ny2)::Array{Float64,2}
  ux = zeros(nx2,ny2)::Array{Float64,2}
  uy = zeros(nx2,ny2)::Array{Float64,2}
  sfd = zeros(nt)::Array{Float64,1}

  dx = 1/dx
  dy = 1/dy
  oneOverRho = 1/rho
  # Main loop
  # no need to separate in a function it does not help performances
  for a = 2:nt
    #println("iteration: $a")
    # Inject source funtion

    px[isx,isy] = px[isx,isy] + dt*0.5*fs[a]
    py[isx,isy] = py[isx,isy] + dt*0.5*fs[a]
     for j=1:ny
       for i=1:nx
        # Update px
         diffop = (ux[i+1,j] - ux[i,j])*dx
         pmlop = qx[i+1,j+1]*px[i+1,j+1]
         px[i+1,j+1] -= dt*(pmlop + rho*vp[i+1,j+1]*diffop)

        # Update py
         diffop = (uy[i,j+1] - uy[i,j])*dy
         pmlop = qy[i+1,j+1]*py[i+1,j+1]
         py[i+1,j+1] -= dt*(pmlop + rho*vp[i+1,j+1]*diffop)

        # Update ux
         diffop = (px[i+1,j+1] - px[i,j+1] + py[i+1,j+1] - py[i,j+1])*dx
         pmlop = 0.5*(qx[i+1,j+1]+qx[i,j+1])*ux[i,j]
         ux[i,j] -= dt*oneOverRho*(pmlop + diffop)

        # Update uy
         diffop = (px[i+1,j+1] - px[i+1,j] + py[i+1,j+1] - py[i+1,j])*dy
         pmlop = 0.5*(qy[i+1,j+1]+qy[i+1,j])*uy[i,j]
         uy[i,j] -= dt*oneOverRho*(pmlop + diffop)
      end
    end
    sfd[a] = px[irx-1,iry-1] + py[irx-1,iry-1]
  end

  return sfd
end

@time fdtd2ds(0)
Profile.clear()
Profile.init(10^9,0.001)
@profile (for i=1:10; fdtd2ds(0); end;)
Profile.print()
#@time fdtd2ds(0)
#tfd = t*1000
#clf
#using Gadfly
#plot(tfd,rec)
#hold on

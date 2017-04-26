
function make_source(source,fc,nt,dt)
  coeff::Float64=4
  tsour::Float64=1/fc
  t::Array{Float64,1}= collect(0:nt-1)*dt
  lsour::Float64=tsour/dt
  t0::Float64=tsour*1.5
  T0::Float64=1/t0
  tau::Array{Float64,1} = zeros(nt)
  @simd for l=1:nt
    tau[l]=pi*(t[l]-t0)*T0
  end
  #tau=pi*(t-t0)/T0
  fs::Array{Float64,1} = zeros(nt)
  if source==0
    @simd for k=1:nt
      fs[k] =(1-coeff*tau[k]*tau[k])*exp(-2*tau[k]*tau[k])
    end
  else
    fs=source
  end
  return fs
end

function fdtd2ds(source, npml_in=80, pmlfac_in=1.2, pmlexp_in=1.5)
  # fdtd2ds: finitite differences time domain 2ds

  # Define velocity model
  const nx::Int = 400
  const ny::Int = 400
  const nx2::Int = nx+1
  const ny2::Int = ny+1
  const vp0::Int = 2000
  vp::Array{Float64,2} = vp0*ones(nx2,ny2)

  #vp = vp .* vp
  @simd for j=1:ny2
    @simd for i=1:nx2
      @inbounds vp[i,j]=vp[i,j]*vp[i,j]
    end
  end
  const rho::Float64 = 1500
  # Define source position and receiver position
  const isx::Int = 200
  const isy::Int = 200
  const irx::Int = 200
  const iry::Int = 300

  # source and fd parameters fd= finite difference
  const nt::Int = 3000 #better is set at 3000
  const fc::Float64 = 1000
  dx::Float64 = vp0/(30*fc)
  dy::Float64 = dx
  dt::Float64 = dx/vp0*1/(sqrt(2)*2*pi)
  #print("dx : $dx , dt : $dt , nt : $nt")

  # Define source-time function
  fs = make_source(source,fc,nt,dt)

  # Define PML (perfectly matched layer boundaries) properties
  npml::Int = npml_in
  npmll::Float64 = npml
  pmlfac::Float64 = pmlfac_in
  pmlexp::Float64 = pmlexp_in
  qx = zeros(nx2,ny2)::Array{Float64,2}# inutile
  qy = zeros(nx2,ny2)::Array{Float64,2}# inutile
  qx = zeros(nx,ny)::Array{Float64,2}
  qy = zeros(nx,ny)::Array{Float64,2}
  @simd for i=1:ny
    @simd for a = 1:npml
      @inbounds qx[a,i] = pmlfac*(npmll-a)^pmlexp # left explained by the transposition of the matrices when plotted  x is vertical axis in the code and y is horizontal, the trickt to change it is to flip the matrices when plotting
      @inbounds qx[nx-a+1,i] = pmlfac*(npmll-a)^pmlexp # right
    end
  end
  @simd for i=1:nx
    @simd for a=1:npml
      @inbounds qy[i,a] = pmlfac*(npmll-a)^pmlexp # thp
      @inbounds qy[i,ny-a+1] = pmlfac*(npmll-a)^pmlexp # bottom
    end
  end

  qx = hcat(qx[:,1], qx) # pour que ça fasse la taille nx2, ny2
  qx = vcat(qx[1,:], qx)
  qy = hcat(qy[:,1], qy)
  qy = vcat(qy[1,:], qy)



  # Initialize fields
  px = zeros(nx2,ny2)::Array{Float64,2}
  py = zeros(nx2,ny2)::Array{Float64,2}
  ux = zeros(nx2,ny2)::Array{Float64,2}
  uy = zeros(nx2,ny2)::Array{Float64,2}
  sfd1 = zeros(nt)::Array{Float64,1}
  sfd2 = zeros(nt)::Array{Float64,1}

  snap = zeros(nx2,ny2,nt)

  dx = 1/dx
  dy = 1/dy
  oneOverRho::Float64 = 1/rho
  # Main loop
  # no need to separate in a function it does not help performances
  for a = 2:nt
    #println("iteration: $a")
    # Inject source funtion

    px[isx,isy] = px[isx,isy] + dt*0.5*fs[a]
    py[isx,isy] = py[isx,isy] + dt*0.5*fs[a]
    @simd for j=1:ny
      @simd for i=1:nx
        # Update px
        @inbounds diffop = (ux[i+1,j] - ux[i,j])*dx
        @inbounds pmlop = qx[i+1,j+1]*px[i+1,j+1]
        @inbounds px[i+1,j+1] = px[i+1,j+1] - dt*(pmlop + rho*vp[i+1,j+1]*diffop)

        # Update py
        @inbounds diffop = (uy[i,j+1] - uy[i,j])*dy
        @inbounds pmlop = qy[i+1,j+1]*py[i+1,j+1]
        @inbounds py[i+1,j+1] = py[i+1,j+1] - dt*(pmlop + rho*vp[i+1,j+1]*diffop)

        # Update ux
        @inbounds diffop = (px[i+1,j+1] - px[i,j+1] + py[i+1,j+1] - py[i,j+1])*dx
        @inbounds pmlop = 0.5*(qx[i+1,j+1]+qx[i,j+1])*ux[i,j]
        @inbounds ux[i,j] = ux[i,j] - dt*oneOverRho*(pmlop + diffop)

        # Update uy
        @inbounds diffop = (px[i+1,j+1] - px[i+1,j] + py[i+1,j+1] - py[i+1,j])*dy
        @inbounds pmlop = 0.5*(qy[i+1,j+1]+qy[i+1,j])*uy[i,j]
        @inbounds uy[i,j] = uy[i,j] - dt*oneOverRho*(pmlop + diffop)
      end
    end
    #snap[:,:,a] = px + py
    sfd1[a] = px[irx,300] + py[irx,300]
    sfd2[a] = px[irx,100] + py[irx,100]

    #if (a % 100 == 0)
    #  println("a: $a/$nt")
    #  plt.imshow(px + py)
    #  plt.draw()
    #elseif a == 2
    #    println("salut")
    #    plt.ion()
    #    plt.imshow(px+py)
    #    plt.show()
    #end
  end

  return sfd1, sfd2
end


##########
## MAIN ##
##########

using PyCall
#@pyimport matplotlib.pyplot as plt

println("repeating 3x the code so that the \"interpreter has warmed up\"")
trace3, trace4 = @time fdtd2ds(0, 30, 50,2)
@time fdtd2ds(0, 30, 50,2)
@time fdtd2ds(0, 30, 50,2)
#
#traces_right = zeros(7,3000)
#traces_left = zeros(7,3000)
#
#for i=1:7
#    println("nmpl: $(i*10)")
#    traces_right[i,:], traces_left[i,:] = fdtd2ds(0, i*10)
#end

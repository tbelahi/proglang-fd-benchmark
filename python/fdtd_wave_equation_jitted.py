import numpy as np
from numba.decorators import autojit, jit
import time
import matplotlib.pyplot as plt

def fdtd2D():
  # define velocity model
  nx = 400
  ny = 400
  nx2 = nx+1
  ny2 = ny+1
  vp0 = 2000.
  vp = vp0*np.ones((nx,ny))

  vp = vp**2
  #vp[ny/2:,:]=4*vp[ny/2:,:]
  rho = 1500.

  # define source and receiver position
  isx = 200
  isy = 100
  irx = 200
  iry = 300

  # source and FD parameters
  nt = 3000
  fc = 1000.
  dx = vp0/(30.*fc) # 30 points per wavelength
  dy = dx
  dt = dx/(vp0*np.sqrt(2)*2*np.pi) # Courant criterion et un facteur 2*pi

  # define source-time function, ricker wavelet
  tsour = 1/fc
  t = np.arange(nt)*dt
  lsour = tsour/dt
  t0 = tsour*1.5
  T0 = tsour*1.5
  tau = np.pi*(t-t0)/T0
  a=4
  fs = (1-a*tau**2)*np.exp(-2*tau**2)

  # define PML properties
  npml = 30
  pmlfac = 50
  pmlexp = 2
  qx = np.zeros((ny,nx))
  qy = np.zeros((ny,nx))
  for a in range(npml):
    qx[:,a] = pmlfac*(npml-a)**pmlexp # left
    qx[:,nx-1-a] = pmlfac*(npml-a)**pmlexp # right
    qy[a,:] = pmlfac*(npml-a)**pmlexp # top
    qy[ny-1-a,:] = pmlfac*(npml-a)**pmlexp # bottom
  # give qx and qy the dimesion (ny2,nx2)
  qx = np.vstack((qx[0,:],qx))
  qx = np.hstack((qx[:,0].reshape((nx2,1)),qx))
  qy = np.vstack((qy[0,:],qy))
  qy = np.hstack((qy[:,0].reshape((nx2,1)),qy))

  # initialize fields
  px = np.zeros((ny2,nx2))
  py = np.zeros((ny2,nx2))
  ux = np.zeros((ny2,nx2))
  uy = np.zeros((ny2,nx2))
  sfd = np.empty((nt,))

  # main loop: looping over timesteps
  looping(qx,qy,px,py,ux,uy,sfd,fs,nt,nx,ny,nx2,ny2,irx,iry,isx,isy,dx,dy,dt,rho,vp)

  return sfd

@jit
def looping(qx,qy,px,py,ux,uy,sfd,fs,nt,nx,ny,nx2,ny2,irx,iry,isx,isy,dx,dy,dt,rho,vp):
 for a in range(1,nt):

    # inject source function
    px[isy,isx] = px[isy,isx] + dt*0.5*fs[a]
    py[isy,isx] = py[isy,isx] + dt*0.5*fs[a]

    for j in range(0,ny-1):
        for i in range(0,nx-1):
            # Update px
            diffop = (ux[i+1,j] - ux[i,j])/dx
            pmlop = qx[i+1,j+1]*px[i+1,j+1]
            px[i+1,j+1] = px[i+1,j+1] - dt*(pmlop + rho*vp[i+1,j+1]*diffop)

            # Update py
            diffop = (uy[i,j+1] - uy[i,j])/dy
            pmlop = qy[i+1,j+1]*py[i+1,j+1]
            py[i+1,j+1] = py[i+1,j+1] - dt*(pmlop + rho*vp[i+1,j+1]*diffop)

            # Update ux
            diffop = (px[i+1,j+1] - px[i,j+1] + py[i+1,j+1] - py[i,j+1])/dx
            pmlop = 0.5*(qx[i+1,j+1]+qx[i,j+1])*ux[i,j]
            ux[i,j] = ux[i,j] - dt/rho*(pmlop + diffop)

            # Update uy
            diffop = (px[i+1,j+1] - px[i+1,j] + py[i+1,j+1] - py[i+1,j])/dy
            pmlop = 0.5*(qy[i+1,j+1]+qy[i+1,j])*uy[i,j]
            uy[i,j] = uy[i,j] - dt/rho*(pmlop + diffop)

    sfd[a] = px[iry-1,irx-1] + py[iry-1,irx-1]

    #if a%100 == 0:
    #  plt.imshow(px+py)
    #  plt.show()


def main():
  t0 = time.time()
  sfd = fdtd2D()
  print 'ellapsed time to run with (brutally used) numba; 1st pass: ' + str(time.time()-t0) + " seconds"
  t0 = time.time()
  sfd = fdtd2D()
  print 'ellapsed time to run with (brutally used) numba; 2nd pass: ' + str(time.time()-t0) + " seconds"
  t0 = time.time()
  sfd = fdtd2D()
  print 'ellapsed time to run with (brutally used) numba; 3rd pass: ' + str(time.time()-t0) + " seconds"


if __name__=='__main__':
  main()

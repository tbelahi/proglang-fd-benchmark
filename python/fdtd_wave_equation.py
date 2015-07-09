#!/usr/bin/env python -tt
#-*- coding: utf-8 -*-

"""" Finite difference time-domain solver for the acoustic wave equation, based on a staggered
grid approach, using Euler's equation and Continuity equation"""

#**********************************
# Thomas Belahi, IPGP, Oct 2013
#
# based on Hansruedi Maurer .m
# file, he provided for modeling
# course at ETH in 2011
# staggered grid approach to acous-
# tic equation
#
# convention:
# - x: horizontal
# - y: vertical, pointing down
#**********************************

import numpy as np
#import matplotlib.pyplot as plt
import time
from numba import jit

@jit
def fdtd2D():
  # define velocity model
  nx = 400
  ny = 400
  nx2 = nx+1
  ny2 = ny+1
  vp0 = 2000.
  vp = vp0*np.ones((nx,ny))

  vp = vp**2
  rho = 1500.

  # define source and receiver position
  isx = 200
  isy = 100
  irx = 200
  iry = 300

  # source and FD parameters
  nt = 3000
  #nt = 300
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
  for a in range(1,nt):

    # inject source function
    px[isy,isx] = px[isy,isx] + dt*0.5*fs[a]
    py[isy,isx] = py[isy,isx] + dt*0.5*fs[a]

    # update px
    diffop = (ux[:ny,1:nx2]-ux[:ny,0:nx])/dx
    pmlop = qx[1:ny2,1:nx2]*px[1:ny2,1:nx2]
    px[1:ny2,1:nx2] = px[1:ny2,1:nx2] - dt*(pmlop+rho*vp*diffop)

    # update py
    diffop = (uy[1:ny2,:nx]-uy[:ny,:nx])/dy
    pmlop = qy[1:ny2,1:nx2]*py[1:ny2,1:nx2]
    py[1:ny2,1:nx2] = py[1:ny2,1:nx2] - dt*(pmlop+rho*vp*diffop)

    # Update ux
    diffop = (px[1:ny2,1:nx2] - px[1:ny2,:nx] + py[1:ny2,1:nx2] - py[1:ny2,:nx])/dx
    pmlop = 0.5*(qx[1:ny2,1:nx2]+qx[1:ny2,:nx])*ux[:nx,:ny]
    ux[:nx,:ny] = ux[:nx,:ny] - dt/rho*(pmlop + diffop)

    #update uy
    diffop = (px[1:ny2,1:nx2] - px[:ny,1:nx2] + py[1:ny2,1:nx2] - py[:ny,1:nx2])/dy
    pmlop = 0.5*(qy[1:ny2,1:nx2]+qy[:ny,1:nx2])*uy[:ny,:nx]
    uy[:nx,:ny] = uy[:nx,:ny] - dt/rho*(pmlop + diffop)

    sfd[a] = px[iry-1,irx-1] + py[iry-1,irx-1]

  return sfd



def main():
  t0 = time.time()
  sfd = fdtd2D()
  print 'ellapsed time to run wit (brutally used) numba; 1st pass: ' + str(time.time()-t0) + " seconds"
  t0 = time.time()
  sfd = fdtd2D()
  print 'ellapsed time to run wit (brutally used) numba; 2nd pass: ' + str(time.time()-t0) + " seconds"
  t0 = time.time()
  sfd = fdtd2D()
  print 'ellapsed time to run wit (brutally used) numba; 3nd pass: ' + str(time.time()-t0) + " seconds"
  #plt.plot(sfd)
  #plt.show()

if __name__=='__main__':
  main()

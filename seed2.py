
import sys;
import numpy as np;
import tfmodule as tf;
from scipy.special import factorial as fac;
from scipy.integrate import simps;
from pathlib import Path;
from math import sqrt;

from numba import jit, prange, int32, uint32, uint64, int64, float64, complex128;




"""
=============================================================================


    GENERATE INITIAL CONDITION DATA
    -------------------------------

    Implementation of useful functions to generate initial conditions for
    the main program, to solve the time dependent GP equation



    CALL IN COMMAND LINE
    --------------------

    python seed.py xi xf Nx yi yf Ny Seed_Name extra_params

    The data is generated in the grid domain [xi,xf] x [yi,yf] with Nx
    discretized points along  x-direction and  Ny along y-direction. A
    'Seed_Name' is required and must be one of the following functions
    in this script. Finally, the files are created in input folder and
    the user must open and change the '*_eq.dat' file to change the eq
    parameters as desired


=============================================================================
"""





def VortexLattice(x,y,den,n):
    n = int(n);
    S = np.zeros([y.size,x.size],dtype=np.complex128);
    S[:,:] = np.sqrt(den[:,:]);

    # phase-noise among the modes of ang. mom. added
    noise = np.pi * (np.random.random(int(n)) - 0.5) / 0.5;
    # grid-noise for every point in domain
    gnoise = (np.random.random(S.shape)-0.5) * 0.2;
    # weight of modes
    w = (1.0*np.ones(n) / fac(np.arange(1,n+1)))**(0.58);

    # Compiled function to sum up all modes
    VortexLatticeAux(x.size,y.size,x,y,n,noise,gnoise,w,den,S);

    # additional grid-random phase
    S = S * np.exp(1.0j*(np.random.random(S.shape)-0.5)*np.pi/4);

    # normalize to 1
    abs2 = abs(S)**2;
    intx = np.zeros(y.size,dtype=np.float64);
    for i in range(y.size): intx[i] = simps(abs2[i], dx = x[1]-x[0]);

    return S / sqrt(simps(intx, dx = y[1]-y[0]));





@jit((uint32,uint32,float64[:],float64[:],uint32,float64[:],float64[:,:],
      float64[:],float64[:,:],complex128[:,:]),
      nopython=True, nogil=True)
def VortexLatticeAux(nx,ny,x,y,n,noise,gnoise,w,den,S):
    Lx = x[-1] - x[0];
    Ly = y[-1] - y[0];
    sigx = 0.33 * sqrt(Lx);
    sigy = 0.33 * sqrt(Ly);
    midx = (x[-1] + x[0]) / 2;
    midy = (y[-1] + y[0]) / 2;
    expargx = 0;
    expargy = 0;
    ph = 0 + 0.0j;
    for i in prange(n):
        for k in prange(ny):
            expargy = - ((y[k] - midy) / sigy)**2
            for l in prange(nx):
                expargx = - ((x[l] - midx) / sigx)**2
                ph = np.exp(1.0j * noise[i] + expargx + expargy)
                # Ang. mom. mode
                rl = (x[l] + 1.0j*y[k])**(i)
                newmode = rl*(1+gnoise[k,l])*w[i]*sqrt(den[k,l])*ph
                S[k,l] = S[k,l] + newmode





# READ DOMAIN GRID INFORMATION
x1 = np.float64(sys.argv[1]);
x2 = np.float64(sys.argv[2]);
nx = int(sys.argv[3]);
y1 = np.float64(sys.argv[4]);
y2 = np.float64(sys.argv[5]);
ny = int(sys.argv[6]);

# EXTRA PARAMETERS TO DEFINE BACKGROUND TF PROFILE
params = []
for i in range(7, len(sys.argv)): params.append(float(sys.argv[i]))
params = tuple(params)

# CHOOSE TRAP PARAMETERS AND INTERACTION FOR THOMAS FERMI 
Ng = params[0]
rot = params[1]
a = params[2]
w = params[3]
b = params[4]

# GENERATE GRID POINTS AND TRAP AT GRID POINTS
x = np.linspace(x1,x2,nx)
y = np.linspace(y1,y2,ny)
V = tf.quarticXQuadraticXY(x,y,w,a,b)

# CHEMICAL POTENTIAL AND FINALLY TF DENSITY PROFILE
mu = tf.searchNG(x,y,Ng,rot,V)
den = tf.den(x,y,mu,rot,V) / Ng

S = VortexLattice(x,y,den,100).reshape(nx*ny);

seedName = 'rotating'

##  Record initial trial function
folder = str(Path.home()) + '/programs/2d-bec/input/';
np.savetxt(folder + seedName + '_init.dat', S.T, fmt='%.15E');

##  Configure domain. Additional default information is setup
##  for time, the step size and quantity, that is  recomended
##  to be adjusted as needed
f = open(folder + seedName + '_domain.dat', "w");
f.write("%.10f %.10f %d %.10f %.10f %d" % (x1, x2, nx, y1, y2, ny));
f.write(" 0.0002 50000");
f.close();

##  Default values are set for equation parameters. One must change
##  for the desired problem
f = open(folder + seedName + '_eq.dat', "w");
if (seedName == 'stationary') : f.write("-0.5 0.0 10.0 1.0 1.0 0.0 0.0");
else : f.write("-0.5 %.14lf %.2lf 1.0 1.0 0.0 0.0 0.0" % (rot,Ng));
f.close();

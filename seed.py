
import sys;
import numpy as np;
from scipy.special import factorial as fac;
from scipy.integrate import simps;
from pathlib import Path;
from math import sqrt;





"""
=============================================================================


    GENERATE INITIAL CONDITION DATA
    -------------------------------

    Given a keyword generate one of some pre-defined initial condition
    data, like solitons, or others yet to be implemented. This initial
    condition can be take to be time evolved or as initial guess  to a
    method to obtain stationary solution for Gross-Pitaevskii equation



    CALL IN COMMAND LINE
    --------------------

    python generate_init.py x1 x2 M Id extra_params

    where [x1, x2] is the position domain and M the number of
    discretized intervals. The number of points therefore  is
    M + 1 in an vector holding each position value. Id is the
    function Identification to generate data and extra_params
    to be able to evalute it.


=============================================================================
"""



def Random(x, y):

    S = np.zeros([y.size,x.size],dtype=np.complex128);

    real = 2 * (np.random.random(S.shape) - 0.5);
    imag = 2 * (np.random.random(S.shape) - 0.5);

    # Localized by a Gaussian like-shape
    sigx = 0.15 * (x[-1] - x[0]);
    sigy = 0.15 * (y[-1] - y[0]);
    midx = (x[-1] + x[0]) / 2;
    midy = (y[-1] + y[0]) / 2;

    expargx = - ((x - midx) / sigx)**2;
    expargy = - ((y - midy) / sigy)**2;

    for k in range(y.size):

        for l in range(x.size):
            ph = expargx[l] + expargy[k];
            r = (x[l] + 1.0j * y[k])
            S[k,l] = (0.6 + r) * np.exp(ph) * (real[k,l] + 1.0j * imag[k,l]);

    abs2 = abs(S)**2;
    intx = np.zeros(y.size,dtype=np.float64);
    for i in range(y.size): intx[i] = simps(abs2[i], dx = x[1]-x[0]);

    return S / sqrt(simps(intx, dx = y[1]-y[0]));





def VortexLattice(x, y, n):
    n = int(n);
    Lx = x[-1] - x[0];
    Ly = y[-1] - y[0];

    S = np.zeros([y.size,x.size],dtype=np.complex128);

    # generate random numbers in the range [-1,1]
    noise = np.pi * (np.random.random(int(n)) - 0.5) / 0.5;
    # generate noise point-wise in the grid
    gnoise = (np.random.random(S.shape)-0.5) * 0.2;

    # Localized by a Gaussian like-shape
    sigx = 0.33 * sqrt(Lx);
    sigy = 0.33 * sqrt(Ly);
    midx = (x[-1] + x[0]) / 2;
    midy = (y[-1] + y[0]) / 2;

    w = np.sqrt(1.0*np.ones(n) / fac(np.arange(1,n+1)));
    expargx = - ((x - midx) / sigx)**2;
    expargy = - ((y - midy) / sigy)**2;

    for i in range(int(n)):

        for k in range(y.size):

            for l in range(x.size):
                ph = 1.0j * noise[i] + expargx[l] + expargy[k];
                rl = (x[l] + 1.0j*y[k])**(i+1);

                S[k,l] = S[k,l] + (rl + gnoise[k,l])*w[i]*np.exp(ph);

    S = S * np.exp(1.0j*(np.random.random(S.shape)-0.5)*np.pi/4);

    abs2 = abs(S)**2;
    intx = np.zeros(y.size,dtype=np.float64);
    for i in range(y.size): intx[i] = simps(abs2[i], dx = x[1]-x[0]);

    return S / sqrt(simps(intx, dx = y[1]-y[0]));




def stationary(x, y):

    S = np.zeros([y.size,x.size],dtype=np.complex128);

    # Localized by a Gaussian like-shape
    sigx = 0.14 * (x[-1] - x[0]);
    sigy = 0.14 * (y[-1] - y[0]);
    midx = (x[-1] + x[0]) / 2;
    midy = (y[-1] + y[0]) / 2;

    expargx = - ((x - midx) / sigx)**2;
    expargy = - ((y - midy) / sigy)**2;

    noise = np.pi * (np.random.random([y.size,x.size]) - 0.5) / 0.5;

    for k in range(y.size):
        for l in range(x.size):
            ph = 1.0j * noise[k,l] + expargx[l] + expargy[k];
            S[k,l] = np.exp(ph);

    abs2 = abs(S)**2;
    intx = np.zeros(y.size,dtype=np.float64);
    for i in range(y.size): intx[i] = simps(abs2[i], dx = x[1]-x[0]);

    return S / sqrt(simps(intx, dx = y[1]-y[0]));





#              Domain discretization parameters               #
#              --------------------------------               #





x1 = np.float64(sys.argv[1]);
x2 = np.float64(sys.argv[2]);
nx = int(sys.argv[3]);

y1 = np.float64(sys.argv[4]);
y2 = np.float64(sys.argv[5]);
ny = int(sys.argv[6]);

seedName = sys.argv[7];

x = np.linspace(x1, x2, nx);
y = np.linspace(y1, y2, ny);

if (seedName == 'rotating') : S = VortexLattice(x,y,35).reshape(nx*ny);
elif (seedName == 'stationary') : S = NoRotation(x,y).reshape(nx*ny);
elif (seedName == 'random') : S = Random(x,y).reshape(nx*ny);
else : print('\nSeed name %s not recognized\n\n' % seedName);

folder = str(Path.home()) + '/programs/2d-bec/input/';

np.savetxt(folder + seedName + '_init.dat', S.T, fmt='%.15E');

f = open(folder + seedName + '_domain.dat', "w");

f.write("%.10f %.10f %d %.10f %.10f %d" % (x1, x2, nx, y1, y2, ny));
f.write(" 0.0005 50000"); # trial values for time step size and quantity

f.close();

f = open(folder + seedName + '_eq.dat', "w");

if (seedName == 'stationary') : f.write("-0.5 0.0 10.0 1.0 1.0 0.0 0.0");
else : f.write("-0.5 0.8 50.0 1.0 1.0 0.0 0.0");

f.close();

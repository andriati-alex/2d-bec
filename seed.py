
import sys;
import numpy as np;
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





def VortexLattice(x, y, n):
    n = int(n);

    S = np.zeros([x.size,y.size],dtype=np.complex128);

    # generate random numbers in the range [-1,1]
    # noisex = (np.random.random(int(n)) - 0.5) / 0.5;
    # noisey = (np.random.random(int(n)) - 0.5) / 0.5;

    # Localized by a Gaussian like-shape
    sigx = 0.16 * (x[-1] - x[0]);
    sigy = 0.16 * (y[-1] - y[0]);
    midx = (x[-1] + x[0]) / 2;
    midy = (y[-1] + y[0]) / 2;

    # Momenta with some noise
    # kx = (np.arange(0, int(n)) * noisex) * 2 * np.pi / sigx;
    # ky = (np.arange(0, int(n)) * noisey) * 2 * np.pi / sigy;

    w = (1.0 / np.arange(1, n + 1)**1.5);
    expargx = - ((x - midx) / sigx)**2;
    expargy = - ((y - midy) / sigy)**2;

    for i in range(int(n)):

        noise = np.pi * (np.random.random([x.size,y.size]) - 0.5) / 0.5;
        for k in range(x.size):

            for l in range(y.size):
                ph = 1.0j * noise[k,l] + expargx[k] + expargy[l];
                rl = 1.0;
                for j in range((i + 1)%2): rl = rl * (x[k] + 1.0j*y[l]);

                S[k,l] = S[k,l] + rl*w[i]*np.exp(ph);

    abs2 = abs(S)**2;
    intx = np.zeros(y.size,dtype=np.float64);
    for i in range(y.size):
        intx[i] = simps(abs2[:,i], dx = x[1]-x[0]);

    return S / sqrt(simps(intx, dx = y[1]-y[0]));





def NoRotation(x, y):

    S = np.zeros([x.size,y.size],dtype=np.complex128);

    # generate random numbers in the range [-1,1]
    # noisex = (np.random.random(int(n)) - 0.5) / 0.5;
    # noisey = (np.random.random(int(n)) - 0.5) / 0.5;

    # Localized by a Gaussian like-shape
    sigx = 0.14 * (x[-1] - x[0]);
    sigy = 0.14 * (y[-1] - y[0]);
    midx = (x[-1] + x[0]) / 2;
    midy = (y[-1] + y[0]) / 2;

    # Momenta with some noise
    # kx = (np.arange(0, int(n)) * noisex) * 2 * np.pi / sigx;
    # ky = (np.arange(0, int(n)) * noisey) * 2 * np.pi / sigy;

    expargx = - ((x - midx) / sigx)**2;
    expargy = - ((y - midy) / sigy)**2;

    noise = np.pi * (np.random.random([x.size,y.size]) - 0.5) / 0.5;

    for k in range(x.size):
        for l in range(y.size):
            ph = 1.0j * noise[k,l] + expargx[k] + expargy[l];
            S[k,l] = np.exp(ph);

    abs2 = abs(S)**2;
    intx = np.zeros(y.size,dtype=np.float64);
    for i in range(y.size):
        intx[i] = simps(abs2[:,i], dx = x[1]-x[0]);

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

if (seedName == 'rotating') : S = VortexLattice(x,y,3).reshape(1,nx*ny)[0];
elif (seedName == 'stationary') : S = NoRotation(x,y).reshape(1,nx*ny)[0];
else : print('\nSeed name %s not recognized\n\n' % seedName);

folder = str(Path.home()) + '/AndriatiLibrary/2d-bec/input/';

np.savetxt(folder + seedName + '_init.dat', S.T, fmt='%.15E');

f = open(folder + seedName + '_domain.dat', "w");

f.write("%.10f %.10f %d %.10f %.10f %d" % (x1, x2, nx, y1, y2, ny));
f.write(" 0.0005 50000"); # trial values for time step size and quantity

f.close();

f = open(folder + seedName + '_eq.dat', "w");

if (seedName == 'stationary') :
    f.write("-0.5 0.0 10.0 1.0 1.0 0.0");
else :
    f.write("-0.5 0.8 50.0 1.0 1.0 0.0");

f.close();

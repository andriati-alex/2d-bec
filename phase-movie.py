import sys;
import numpy as np;
import matplotlib.pyplot as plt;

from pathlib import Path;
from matplotlib.animation import FuncAnimation;

fname_prefix = sys.argv[1];
folder = str(Path.home()) + '/programs/2d-bec/output/';
full_path = folder + fname_prefix;

S = np.loadtxt(full_path + '_orb_realtime.dat',dtype=np.complex128);
conf = np.loadtxt(full_path + '_conf_realtime.dat');
t = np.loadtxt(full_path + '_time.dat');

# compute the phase previously to do not compute in frame deneration
ph = np.empty(S.shape);
for i in range(S.shape[0]) : ph[i] = np.arctan2(S[i].imag,S[i].real);

xi = conf[0];
xf = conf[1];
yi = conf[3];
yf = conf[4];
nx = int(conf[2]);
ny = int(conf[5]);
Nsteps = t.size;

animTime = 10000 # in miliseconds = 10 seconds

fig = plt.figure();

im = plt.imshow(ph[0].reshape(ny,nx),interpolation='hanning',cmap='seismic',
        extent=[xi,xf,yi,yf],aspect='equal',origin='lower');

# initialization function: plot the background of each frame
def init():
    im.set_data(ph[0].reshape(ny,nx));
    return [im]

# animation function. This is called sequentially
def animate(i):
    im.set_data(ph[i].reshape(ny,nx));
    return [im]

anim = FuncAnimation(fig, animate, init_func=init,
                     frames=Nsteps, interval=animTime/Nsteps, blit=True);

plt.show();

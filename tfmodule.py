import numpy as np

from math import pi, sqrt
from scipy.integrate import simps
from numba import jit, prange, int32, float64

def den(x,y,mu,rot,pot):
    """ COMPUTE DENSITY ACCORDING TO THOMAS-FERMI """
    nx = x.size
    ny = y.size
    density = np.zeros([ny,nx],dtype=np.float64)
    NUMBA_den(nx,ny,x,y,mu,rot,pot,density)
    return density

@jit((int32,int32,float64[:],float64[:],float64,float64,float64[:,:],
      float64[:,:]),nopython=True,nogil=True)
def NUMBA_den(nx,ny,x,y,mu,rot,pot,density):
    for j in prange(ny):
        for i in prange(nx):
            r = sqrt(y[j]**2 + x[i]**2)
            D = mu - pot[j,i] + 0.5*rot*rot*r*r;
            # Setup value if it is valid according to TF-regime
            if (D > 0): density[j,i] = D

def quarticXQuadraticXY(x,y,w,a,b):
    """ TRAP POTENTIAL """
    xgrid = x.size
    ygrid = y.size
    V = np.zeros([ygrid,xgrid],dtype=np.float64)
    w2 = w*w
    a2 = a*a
    b2 = b*b
    a4 = a*a*a*a
    for j in range(ygrid):
        y2 = y[j]*y[j]
        for i in range(xgrid):
            x2 = x[i]*x[i]
            x4 = x2*x2
            V[j,i] = 0.5*(w2*y2 - a2*x2) + 0.25*b2*x4 + 0.25*a4/b2
    return V

def quarticXYQuadraticXY(x,y,w,a,b):
    """ TRAP POTENTIAL """
    xgrid = x.size
    ygrid = y.size
    V = np.zeros([ygrid,xgrid],dtype=np.float64)
    w2 = w*w
    a2 = a*a
    b2 = b*b
    a4 = a*a*a*a
    for j in range(ygrid):
        y2 = y[j]*y[j]
        for i in range(xgrid):
            x2 = x[i]*x[i]
            r4 = (x2 + y2)*(x2 + y2)
            V[j,i] = 0.5*(w2*y2 - a2*x2) + 0.25*b2*r4 + 0.25*a4/b2
    return V

def searchNG(x,y,Ng,rot,pot,tol=1E-9,guess=35.9764):
    """ SEARCH \mu CORRESPONDING TO THE OTHER FIXED PARAMETERS """
    nx = x.size
    ny = y.size
    dx = (x[nx-1] - x[0]) / (nx - 1)
    dy = (y[ny-1] - y[0]) / (ny - 1)
    # compute at inital arbitrary value
    mu = guess
    density = den(x,y,mu,rot,pot)
    intx = np.zeros(ny,dtype=np.float64)
    for j in range(ny): intx[j] = simps(density[j],dx=dx);
    trial = simps(intx,dx=dy)
    diff = trial - Ng
    # DEFINE THE LIMITS WHERE SEARCH FOR \mu
    # If it is above desired decrease in fixed steps to find minimum
    while (diff > 0):
        mu = mu - 1
        density = den(x,y,mu,rot,pot)
        for j in range(ny): intx[j] = simps(density[j],dx=dx);
        trial = simps(intx,dx=dy)
        diff = trial - Ng
    mumin = mu - 1
    # If it is below desired increase in fixed steps to find maximum
    while (diff < 0):
        mu = mu + 1
        density = den(x,y,mu,rot,pot)
        for j in range(ny): intx[j] = simps(density[j],dx=dx);
        trial = simps(intx,dx=dy)
        diff = trial - Ng
    mumax = mu + 1
    while (abs(diff) > tol):
        if (diff > 0):
            mumax = mu
            mu = mu - (mumax - mumin) / 2
        else:
            mumin = mu
            mu = mu + (mumax - mumin) / 2
        density = den(x,y,mu,rot,pot)
        for j in range(ny): intx[j] = simps(density[j],dx=dx);
        trial = simps(intx,dx=dy)
        diff = trial - Ng
    print("mu = %.7lf | err = %.1E" % (mu,abs(diff)))
    return mu

def critRot(x,y,Ng,a,b,w,pot,tol=1E-6,guess=1.0):
    """ COMPUTE CRITICAL ROTATION TO BREAK IN TWO PIECES """
    xmin = a / b
    rot = guess
    rotmin = 0.0 # no rotation at all
    rotmax = w   # if exceed it the condensate flies appart
    mu = searchNG(x,y,Ng,rot,pot)
    # diff must be zero to fulfill breakup (see the notes)
    diff = 4 * mu / b**2 - xmin**4
    while (abs(diff) > 0):
        print("Rotation : %.4lf" % rot)
        if (diff > 0):
            rotmin = rot
            rot = rot + (rotmax - rotmin) / 2
        else:
            rotmax = rot
            rot = rot - (rotmax - rotmin) / 2
        mu = searchNG(x,y,Ng,rot,pot,guess=1.2*mu)
        diff = 4 * mu / b**2 - xmin**4
        if (rotmax - rotmin < tol): break
        if (rotmin > 0.999*w):
            print("\nCritical rotation exceeded harmonic confinement.\n")
            return w
    print("\nCritical rotation : %.4lf\n" % rot)
    return rot

def critRotQuarticR(x,y,Ng,a,b,w,pot,tol=1E-6,guess=0.5):
    """ COMPUTE CRITICAL ROTATION TO BREAK IN TWO PIECES """
    xmin = a / b
    rot = guess
    rotmin = 0.0 # no rotation at all
    rotmax = 3*w # upper bound for rotation
    mu = searchNG(x,y,Ng,rot,pot)
    # diff must be zero to fulfill breakup (see the notes)
    diff = 4 * mu / b**2 - xmin**4
    while (abs(diff) > 0):
        print("Rotation : %.4lf" % rot)
        if (diff > 0):
            rotmin = rot
            rot = rot + (rotmax - rotmin) / 2
        else:
            rotmax = rot
            rot = rot - (rotmax - rotmin) / 2
        mu = searchNG(x,y,Ng,rot,pot,guess=1.2*mu)
        diff = 4 * mu / b**2 - xmin**4
        if (rotmax - rotmin < tol): break
    print("\nStep Critical rotation : %.4lf\n" % rot)
    rot_step = rot
    if (rot > w):
        rotmin = rot*0.85
        print("\nCenter Hole phase. Increase rotation to breakup\n")
        while((rot**2 - w**2)**2/b**4 > abs(diff)):
            rot = rot * 1.2 # increase 20%
            mu = searchNG(x,y,Ng,rot,pot,guess=mu)
            diff = 4 * mu / b**2 - xmin**4
        # End of while found maximum rotation
        rotmax = rot*1.1
        rot = 0.5*(rotmax + rotmin)
        mu = searchNG(x,y,Ng,rot,pot,guess=mu)
        diff = (rot**2 - w**2)**2/b**4 + (4 * mu / b**2 - xmin**4)
        while (abs(diff) > 0):
            print("Rotation : %.4lf" % rot)
            if (diff > 0):
                rotmin = rot
                rot = rot + (rotmax - rotmin) / 2
            else:
                rotmax = rot
                rot = rot - (rotmax - rotmin) / 2
            mu = searchNG(x,y,Ng,rot,pot,guess=1.2*mu)
            diff = (rot**2 - w**2)**2/b**4 + (4 * mu / b**2 - xmin**4)
            if (rotmax - rotmin < tol): break

    print("\nFinal Critical rotation : %.4lf\n" % rot)
    return rot_step, rot

def critRot_b(x,y,Ng,a,w,binit=0.1,bfinal=1.0,bstep=0.01):
    """ COMPUTE CRITICAL ROTATION FREQ. VARYING QUARTIC PARAMETER """
    bvals = np.arange(bfinal,binit-bstep/2,-bstep)
    critical = np.zeros(bvals.size)
    initGuess = 1.5
    k = 0
    for b in bvals:
        # update potential because a parameter has changed
        V = quarticXQuadraticXY(x,y,w,a,b)
        # find critical rotation frequency
        critical[k] = critRot(x,y,Ng,a,b,w,V,guess=initGuess)
        initGuess = critical[k] * 0.9
        k = k + 1
    return bvals[::-1], critical[::-1]

def critRot_b_QuarticR(x,y,Ng,a,w,binit=0.2,bfinal=1.0,bstep=0.2):
    """ COMPUTE CRITICAL ROTATION FREQ. VARYING QUARTIC PARAMETER """
    bvals = np.arange(bfinal,binit-bstep/2,-bstep)
    critical = np.zeros(bvals.size)
    hole = np.zeros(bvals.size)
    initGuess = w/2
    k = 0
    for b in bvals:
        print("\n\n\nDOING b = %.3lf\n\n\n" % b)
        # update potential because a parameter has changed
        V = quarticXYQuadraticXY(x,y,w,a,b)
        # find critical rotation frequency
        hole[k], critical[k] = critRotQuarticR(x,y,Ng,a,b,w,V,guess=initGuess)
        initGuess = hole[k] * 0.9
        k = k + 1
    return bvals[::-1], hole[::-1], critical[::-1]

def radiusPlus(rot,w,mu,a,b):
    xmin = a / b
    gam = 4 * mu / b**2 - xmin**4
    # avoid numerical computation exactly at odd multiple
    theta1 = 2 * np.pi * np.arange(0,0.25,0.001)
    theta2 = 2 * np.pi * np.arange(0.251,0.75,0.001)
    theta3 = 2 * np.pi * np.arange(0.751,1.0,0.001)
    theta = np.concatenate([theta1,theta2,theta3])
    r = np.zeros(theta.size)
    w2b2 = w**2 / b**2
    rot2b2 = rot**2 / b**2
    xmin2 = xmin**2
    forb = 0
    # Run over the angle values and ignore those which would
    # give complex radius to at the end exclude them
    for i in range(theta.size):
        t = np.tan(theta[i])**2
        inRoot = (xmin2 - w2b2*t + rot2b2*(1+t))**2 + gam
        # Check if the angle there is a real solution for 'r'
        if (inRoot >= 0):
            outRoot = xmin2 - w2b2*t + rot2b2*(1+t)
            if (outRoot + np.sqrt(inRoot) >= 0):
                r[i] = np.sqrt((1 + t) * (outRoot + np.sqrt(inRoot)))
            else:
                forb = forb + 1
                r[i] = -1
        else:
            forb = forb + 1
            r[i] = -1
    # Exclude the forbidden region where there is no real solution
    newr = np.zeros(theta.size - forb)
    newtheta = np.zeros(theta.size - forb)
    k = 0
    for i in range(theta.size):
        if (r[i] >= 0):
            newtheta[k] = theta[i]
            newr[k] = r[i]
            k = k + 1
    return newtheta, newr

def radiusMinus(rot,w,mu,a,b):
    xmin = a / b
    gam = 4 * mu / b**2 - xmin**4
    # avoid numerical computation exactly at odd multiple
    theta1 = 2 * np.pi * np.arange(0,0.25,0.001)
    theta2 = 2 * np.pi * np.arange(0.251,0.75,0.001)
    theta3 = 2 * np.pi * np.arange(0.751,1.0,0.001)
    theta = np.concatenate([theta1,theta2,theta3])
    r = np.zeros(theta.size)
    w2b2 = w**2 / b**2
    rot2b2 = rot**2 / b**2
    xmin2 = xmin**2
    forb = 0
    # Run over the angle values and ignore those which would
    # give complex radius to at the end exclude them
    for i in range(theta.size):
        t = np.tan(theta[i])**2
        inRoot = (xmin2 - w2b2*t + rot2b2*(1+t))**2 + gam
        # Check if the angle there is a real solution for 'r'
        if (inRoot >= 0):
            outRoot = xmin2 - w2b2*t + rot2b2*(1+t)
            if (outRoot - np.sqrt(inRoot) >= 0):
                r[i] = np.sqrt((1 + t) * (outRoot - np.sqrt(inRoot)))
            else:
                forb = forb + 1
                r[i] = -1
        else:
            forb = forb + 1
            r[i] = -1
    # Exclude the forbidden region where there is no real solution
    newr = np.zeros(theta.size - forb)
    newtheta = np.zeros(theta.size - forb)
    k = 0
    for i in range(theta.size):
        if (r[i] >= 0):
            newtheta[k] = theta[i]
            newr[k] = r[i]
            k = k + 1
    return newtheta, newr

def FeynmanVor(x,y,mu,rot,pot):
    """ COMPUTE DENSITY ACCORDING TO THOMAS-FERMI """
    nx = x.size
    ny = y.size
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    hasCloud = np.zeros([ny,nx])
    for j in prange(ny):
        for i in prange(nx):
            r = sqrt(y[j]**2 + x[i]**2)
            D = mu - pot[j,i] + 0.5*rot*rot*r*r;
            # Setup value if it is valid according to TF-regime
            if (D > 0): hasCloud[j,i] = 1
    # compute area
    intx = np.zeros(ny)
    for i in range(ny): intx[i] = simps(hasCloud[i],dx=dx)
    area = simps(intx,dx=dy)
    return area * rot / np.pi

import numpy as np
import matplotlib.pyplot as plt

def intitial_conditions(x):
    
    return np.hstack((np.zeros(17), np.linspace(0.0, 1.0, 7), np.zeros(17)))
    #return np.hstack((np.zeros(17), np.random.rand(7), np.zeros(17)))
    #return np.where(x%1. < 0.5, np.power(np.sin(2*x*np.pi), 2), 0)

def main():

    nx = 40
    nt = 10
    ## courant number is equal to u*deltaT/deltaX ##
    x = np.linspace(0.0, 1.0, nx+1)
    ## set quantities ##
    u = 5.0    
    c = 0.2
    dx = 1.0/nx
    dt = c*dx/u
    t = nt*dt

    ## initial conditions ##
    x = np.linspace(0.0, 1.0, nx+1)
    phi_initial = intitial_conditions(x)

    phi = phi_initial.copy()
    phinew = phi_initial.copy()
    phiold = phi_initial.copy()

    ''' FTFS '''

    ## reset initial conditions ##
    phi = phi_initial.copy()
    phinew = phi_initial.copy()
    phiold = phi_initial.copy()

    for time in range(0,nt):
        for space in range(0,nx,1):
            phinew[space] = phiold[space] - c*(phiold[space+1] - phiold[space])
        # apply periodic boundaries #
        phinew[0] = phiold[0] - c*(phiold[1] - phiold[0])
        phinew[nx] = phinew[0]
        phiold = phi.copy()
        phi = phinew.copy()

    phiFTFS = phi

    ''' FTBS '''

    ## reset initial conditions ##
    phi = phi_initial.copy()
    phinew = phi_initial.copy()
    phiold = phi_initial.copy()

    for time in range(0,nt):
        for space in range(nx,0,-1):
            phinew[space] = phiold[space] - c*(phiold[space] - phiold[space-1])
        # apply periodic boundaries #
        phinew[nx] = phiold[nx] - c*(phiold[nx] - phiold[0])
        phinew[0] = phinew[nx] 
        phiold = phi.copy()
        phi = phinew.copy()

    phiFTBS = phi

    
    ''' CTCS '''

    ## reset initial conditions ##
    phi = phi_initial.copy()
    phinew = phi_initial.copy()
    phiold = phi_initial.copy()

    ## FTCS for the first time-step, looping over space ##
    for space in range(1,nx):
        phi[space] = phiold[space] - 0.5*c*(phi[space+1] - phi[space-1])
    # apply periodic boundaries #
    phi[0] = phiold[0] - 0.5*c*(phiold[1] - phiold[nx-1])
    phi[nx] = phi[0]

    ## do 50, 100, 150 ##
    ## CTCS blows up when it crosses boundaries ##
    ## CTCS for the first time-step, looping over space ##
    for time in range(1,nt):
        for space in range(1,nx):
            phinew[space] = phiold[space] - c*(phi[space+1] - phi[space-1])
        # apply periodic boundaries #
        phinew[0] = phiold[0] - c*(phi[1] - phi[space-1])
        phinew[nx] = phinew[0]
        phiold = phi.copy()
        phi = phinew.copy()

    phiCTCS = phi
    
    ## plot ##
    plt.figure(figsize=(12, 12))
    plt.plot(x, intitial_conditions(x), color ='k', linestyle = ':', label = 'Initial Conditions')
    plt.plot(x, intitial_conditions(x - u*t), color = 'k', label = 'Analytical')
    plt.plot(x, phiCTCS, color = 'b', label = 'CTCS (FTCS first step)')
    #plt.plot(x, phiFTFS, color = 'r', label = 'FTFS')
    #plt.plot(x, phiFTFS, color = 'g', label = 'FTBS', linestyle = '--')
    plt.xlim(xmin = 0.0, xmax = 1.0)
    plt.ylim(ymin = 0.0, ymax = 1.1)
    plt.legend(loc = 'best')
    plt.ylabel('$\phi$')
    plt.title('Linear Advection Schemes after {} Timesteps'.format(nt))
    #plt.axhline(0, linestyle = ':', color = 'k')
    plt.show()
     
main()

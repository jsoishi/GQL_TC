"""
Usage:
  TC_stability.py  --re=<reynolds> --eta=<eta> [--m=<initial_m>] [--nr=<nr>] [--mu=<mu>] [--ar=<Gamma>] [--re2=<re2>] [--re_range=<re_range> --kz_range=<kz_range>]

Options:
  --nr=<nr>                 number of Radial modes [default: 32]
  --re=<reynolds>           Reynolds number for simulation
  --eta=<eta>               Eta - ratio of R1/R2
  --m=<initial_m>           M1 mode to begin initial conditions [default: 0]
  --mu=<mu>                 mu = Omega2/Omega1 [default: 0]
  --re2=<re2>               Re2 [default: None]
  --ar=<Gamma>              Aspect ratio (height/width) [default: 3]
  --re_range=<re_range>     range of Re for crit finder [default: None]
  --kz_range=<kz_range>     range of kz for crit finder [default: None]

"""
import numpy as np
from dedalus import public as de
from eigentools import Eigenproblem, CriticalFinder
import matplotlib.pyplot as plt
from docopt import docopt
import os
from mpi4py import MPI

import logging
logger = logging.getLogger(__name__.split('.')[-1])


args=docopt(__doc__)
nr = int(args['--nr'])
Re1=float(args['--re'])
if args['--re2'] == 'None':
    Re2 = None
else:
    Re2=float(args['--re2'])
eta=float(args['--eta'])
mu = float(args['--mu'])
Gamma = float(args['--ar'])
m = int(args['--m'])
if args['--re_range'] == 'None':
    re_range = (350,450)
else:
    re_range = [float(i) for i in args['--re_range'].split(',')]
if args['--kz_range'] == 'None':
    kz_range = (350,450)
else:
    kz_range = [float(i) for i in args['--kz_range'].split(',')]
plt.style.use('prl')

alpha = 2*np.pi/Gamma
if Re2:
    root_name = "TC_stability_eta_{:e}_Gamma_{:e}_m_{:d}_Re2_{:e}_nr_{:d}".format(eta,Gamma,m,Re2,nr)
    if MPI.COMM_WORLD.rank == 0:
        logger.info("Computing Stability for eta={:e}, m={:d}, Gamma={:e}, Re2 = {:e}, nr = {:d}".format(eta, m, Gamma, Re2, nr))
        logger.info("Computing Stability with  {:f} <= Re_1 <= {:f}, {:f} <= kz <= {:f}".format(re_range[0], re_range[1], kz_range[0], kz_range[1]))
else:
    root_name = "TC_stability__eta_{:e}_Gamma_{:e}_m_{:d}_mu_{:e}_nr_{:d}".format(eta,Gamma,m,mu,nr)
    if MPI.COMM_WORLD.rank == 0:
        logger.info("Computing Stability for eta={:e}, m={:d}, Gamma={:e}, mu = {:e}, nr = {:d}".format(eta, m, Gamma, mu, nr))
"""
delta = R2 - R1
mu = Omega2/Omega1
eta = R1/R2

scale [L] = delta
scale [T] = delta/(R1 Omega1)
scale [V] = R1 Omega1
"""
#derived parameters
R1 = eta/(1. - eta)
R2 = 1./(1-eta)
nu = 1./Re1

variables = ['u','ur','v','vr','w','wr','p']

#domain
r_basis = de.Chebyshev('r', nr, interval=[R1, R2])

bases = [r_basis]
domain = de.Domain(bases, comm=MPI.COMM_SELF) 

#problem
problem = de.EVP(domain, eigenvalue='sigma', variables=variables)

#params into equations
problem.parameters['eta']=eta
if Re2:
    problem.parameters['Re2']=Re2
else:
    problem.parameters['mu'] = mu
problem.parameters['Re1']=Re1
problem.parameters['kz'] = alpha
problem.parameters['m'] = m

#Substitutions

"""
this implements the cylindrical del operators. 
NB: ASSUMES THE EQUATION SET IS PREMULTIPLIED BY A POWER OF r (SEE BELOW)!!!

Lap_s --> scalar laplacian
Lap_r --> r component of vector laplacian
Lap_t --> theta component of vector laplacian
Lap_z --> z component of vector laplacian

"""
problem.substitutions['nu'] = '1/Re1'
if Re2:
    problem.substitutions['mu'] = 'eta*Re2/Re1'
problem.substitutions['A'] = '(1/eta - 1.)*(mu-eta**2)/(1-eta**2)'
problem.substitutions['B'] = 'eta*(1-mu)/((1-eta)*(1-eta**2))'

problem.substitutions['v0'] = 'A*r + B/r'       #background profile? forcing instead of forcing the boundaries
problem.substitutions['dv0dr'] = 'A - B/(r*r)'  #d/dr of background forcing

problem.substitutions['dtheta(f)'] = '1j*m*f'
problem.substitutions['dz(f)'] = '1j*kz*f'
problem.substitutions['dt(f)'] = 'sigma*f'

# assume pre-multiplication by r*r
problem.substitutions['Lap_s(f, f_r)'] = "r*r*dr(f_r) + r*f_r + dtheta(dtheta(f)) + r*r*dz(dz(f))"
problem.substitutions['Lap_r'] = "Lap_s(u, ur) - u - 2*dtheta(v)"
problem.substitutions['Lap_t'] = "Lap_s(v, vr) - v + 2*dtheta(u)"
problem.substitutions['Lap_z'] = "Lap_s(w, wr)"

# momentum equations
problem.add_equation("r*r*dt(u) - nu*Lap_r - 2*r*v0*v + r*v0*dtheta(u) + r*r*dr(p) = 0")
problem.add_equation("r*r*dt(v) - nu*Lap_t + r*r*dv0dr*u + r*v0*u + r*v0*dtheta(v) + r*dtheta(p)  = 0")
problem.add_equation("r*r*dt(w) - nu*Lap_z + r*r*dz(p) + r*v0*dtheta(w) = 0.")

#continuity
problem.add_equation("r*ur + u + dtheta(v) + r*dz(w) = 0")

#Auxillilary equations
problem.add_equation("ur - dr(u) = 0")
problem.add_equation("vr - dr(v) = 0")
problem.add_equation("wr - dr(w) = 0")

#Boundary Conditions
problem.add_bc("left(u) = 0")
problem.add_bc("right(u) = 0")
problem.add_bc("left(v) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")


ep = Eigenproblem(problem)

cf = CriticalFinder(ep, ("kz","Re1"), find_freq=False)
nx = 10#0
ny = 10#0
xpoints = np.linspace(kz_range[0], kz_range[1], nx)
ypoints = np.linspace(re_range[0], re_range[1], ny)

grid_file = 'results/'+root_name+'grid'
if os.path.exists(grid_file+'.h5'):
    if MPI.COMM_WORLD.rank == 0:
        cf.load_grid(grid_file+'.h5')
else:
    cf.grid_generator((xpoints,ypoints))
    cf.save_grid(grid_file)
if MPI.COMM_WORLD.rank == 0:
    try:
        crits = cf.crit_finder(maxiter=300)
        print(crits)
        logger.info("m = {:d}, eta = {:5.4f}: Critical Re = {:5.2f}, kz = {:7.5f}".format(m, eta, crits[1],crits[0]))
    except ValueError:
        
        crits = None
        print("Critical finder failed.")
    pax,cax = cf.plot_crit()
    #pax.collections[0].set_clim(-0.03,-0.08)

    cax.xaxis.set_ticks_position('top')
    cax.xaxis.set_label_position('top')
    pax.contour(cf.parameter_grids[0], cf.parameter_grids[1],cf.evalue_grid.real, levels=(0,), colors='white')
    pax.figure.savefig('../figs/'+root_name+'_growth_rates.png',dpi=300)

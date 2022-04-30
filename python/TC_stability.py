"""
Usage:
  TC_stability.py  --re=<reynolds> --eta=<eta> [--m=<initial_m>] [--nr=<nr>] [--mu=<mu>] [--ar=<Gamma>] 

Options:
  --nr=<nr>        number of Radial modes [default: 32]
  --re=<reynolds>  Reynolds number for simulation
  --eta=<eta>      Eta - ratio of R1/R2
  --m=<initial_m>  M1 mode to begin initial conditions [default: 0]
  --mu=<mu>        mu = Omega2/Omega1 [default: 0]
  --ar=<Gamma>     Aspect ratio (height/width)
"""
import numpy as np
from dedalus import public as de
from eigentools import Eigenproblem, CriticalFinder
import matplotlib.pyplot as plt
from docopt import docopt
import os


args=docopt(__doc__)
nr = int(args['--nr'])
Re1=float(args['--re'])
eta=np.float(args['--eta'])
mu = np.float(args['--mu'])
Gamma = float(args['--ar'])
m = int(args['--m'])

plt.style.use('prl')

alpha = 2*np.pi/Gamma
root_name = "TC_stability_re_{:e}_eta_{:e}_Gamma_{:e}_m_{:d}_mu_{:e}_nr_{:d}".format(Re1,eta,Gamma,m,mu,nr)

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
domain = de.Domain(bases) 

#problem
problem = de.EVP(domain, eigenvalue='sigma', variables=variables)

#params into equations
problem.parameters['eta']=eta
problem.parameters['mu']=mu
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
ep.solve(sparse=False)
print("Re1={:e}; eta={:e}; Gamma={:e}; m={:d}; mu= {:e} nr={:d}".format(Re1,eta,Gamma,m,mu,nr))
print("Max growth rate = {:e}".format(ep.evalues_good.real.max()))


cf = CriticalFinder(ep, ("kz","Re1"), find_freq=False)
nx = 10
ny = 10
xpoints = alpha*np.linspace(1, 2, nx)
ypoints = np.linspace(100,150,ny)

grid_file = 'results/'+root_name+'grid'
if os.path.exists(grid_file+'.h5'):
    cf.load_grid(grid_file+'.h5')
else:
    cf.grid_generator((xpoints,ypoints))
    cf.save_grid(grid_file)
crits = cf.crit_finder()
print("Critical Re = {:5.2f}, kz = {:7.5f}".format(crits[1],crits[0]))
pax,cax = cf.plot_crit()
pax.collections[0].set_clim(0,0.03)
cax.xaxis.set_ticks_position('top')
cax.xaxis.set_label_position('top')
pax.contour(cf.parameter_grids[0], cf.parameter_grids[1],cf.evalue_grid.real, levels=(0,), colors='white')
pax.figure.savefig('../figs/'+root_name+'_growth_rates.png',dpi=300)

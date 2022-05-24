"""
Usage:
  TC_pseudospectrum.py [--re=<reynolds>] [--eta=<eta>] [--m=<initial_m>] [--kz=<kz>] [--mu=<mu>] [--nr=<nr>] 

Options:
  --re=<reynolds>  Reynolds number for simulation [default: 250]
  --eta=<eta>      Eta - ratio of R1/R2 [default: 0.5]
  --m=<initial_m>  theta mode [default: 0]
  --kz=<kz>        z mode [default: 3.14159265359]
  --mu=<mu>        mu = Omega2/Omega1 [default: -1]
  --nr=<nr>        radial modes [default: 32]

delta = R2 - R1
mu = Omega2/Omega1
eta = R1/R2

scale [L] = delta
scale [T] = delta/(R1 Omega1)
scale [V] = R1 Omega1

Default parameters from Hristova et al (2002, Phys. Fluids), scaled to our non-dimensionalization:
Re1 and kz are multiplied by 2.
"""
from docopt import docopt
import numpy as np
from dedalus import public as de
from eigentools import Eigenproblem, CriticalFinder
import matplotlib.pyplot as plt

args=docopt(__doc__)
Re1=float(args['--re'])
eta=float(args['--eta'])
mu = float(args['--mu'])
nr=int(args['--nr'])
m = int(args['--m'])
alpha = float(args['--kz'])

filename = "../figs/TC_pseudospectrum_Re1_{:5.2f}_eta_{:6.4f}_mu_{:6.4f}_nr_{:d}_m_{:d}_kz_{:5.2f}.png".format(Re1, eta, mu, nr, m, alpha)

def energy_norm(Q1, Q2):
    u1 = Q1['u']
    v1 = Q1['v']
    w1 = Q1['w']
    u2 = Q2['u']
    v2 = Q2['v']
    w2 = Q2['w']

    field = (np.conj(u1)*u2 + np.conj(v1)*v2 + np.conj(w1)*w2).evaluate().integrate()
    return field['g'][0]

plt.style.use('prl')

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
problem.parameters['nu']=nu
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

k = 25

psize = 100
ps_real = np.linspace(-4,2,psize)
ps_imag = np.linspace(-3,3,psize)
ep.calc_ps(k, (ps_real, ps_imag), inner_product=energy_norm,maxiter=20)

clevels = np.linspace(-1.8,0,7)
print("contour levels: ", clevels)
plt.figure(figsize=(5,5))
plt.scatter(ep.evalues_good.imag, ep.evalues_good.real)
plt.xlim(-1.5,1.5)
plt.ylim(-2,1)
plt.axhline(0,color='k',alpha=0.3)
plt.contour(ep.ps_imag, ep.ps_real, np.log10(ep.pseudospectrum.T),levels=clevels,colors='k')
plt.xlabel(r"$\omega$")
plt.ylabel(r"$\gamma$")
plt.tight_layout()
plt.savefig(filename, dpi=300)

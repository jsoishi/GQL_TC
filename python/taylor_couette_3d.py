"""
Usage:
  taylor_couette_3d.py --re=<reynolds> --eta=<eta> --m=<initial_m> [--mu=<mu>] [--ar=<Gamma>] [--nr=<nr>] [--nt=<nt>] [--nz=<nz>] [--restart=<restart_file>] --mesh_1=<mesh_1> --mesh_2=<mesh_2> [--GQL=<GQL>] [--Lambda_z=<Lambda_z>] [--Lambda_theta=<Lambda_theta>] [--willis] [--single_mode] [--run_note=<run_note>] [--theta_symmetry] [--perturb_restart] [--stop_time=<stop_time>] taylor_couette_3d.py

Options:
  --re=<reynolds>  Reynolds number for simulation
  --eta=<eta>      Eta - ratio of R1/R2
  --m=<initial_m>  M1 mode to begin initial conditions
  --mu=<mu>        mu = Omega2/Omega1 [default: 0]
  --ar=<Gamma>     Aspect ratio (height/width)
  --nr=<nr>        radial modes [default: 32]
  --nt=<nt>        theta modes [default: 32]
  --nz=<nz>        z modes [default: 32]
  --mesh_1=<mesh_1> First mesh core count
  --mesh_2=<mesh_2> Second mesh core count
  --restart=<restart_file>    Point to a merged snapshots_s1.h5 file [default: None]
  --GQL=<GQL> True or False
  --Lambda_z=<Lambda_z> Specify an integer cutoff to seperate low and high modes for z
  --Lambda_theta=<Lambda_theta> Specify an integer cutoff to seperate low and high modes for theta
  --willis  Use Willis ICs [default: False]  
  --single_mode  Use single mode ICs [default: False]
  --run_note=<run_note>  Note to add to run directory name [default: None]
  --theta_symmetry  Restrict theta to 2pi/m1 [default: False]
  --perturb_restart    perturb on restart [default: False]
  --stop_time=<stop_time>        time to run (will be added to restart time on restart) [default: None]
delta = R2 - R1
mu = Omega2/Omega1
eta = R1/R2

scale [L] = delta
scale [T] = delta/(R1 Omega1)
scale [V] = R1 Omega1

Default parameters from Barenghi (1991, J. Comp. Phys.).
"""

import numpy as np
import h5py
from dedalus.extras import flow_tools
from dedalus import public as de
from dedalus.tools  import post
import time
import logging
from docopt import docopt
import os
import subprocess
from mpi4py import MPI
from GQLProjection import GQLProject
from random_vector_potential import random_vector_potential
from condition import cond_number

logger = logging.getLogger(__name__)
comm=MPI.COMM_WORLD
rank=comm.Get_rank()

args=docopt(__doc__)
Re1=float(args['--re'])
eta=float(args['--eta'])
mu = float(args['--mu'])
Gamma = float(args['--ar'])
GQL = args['--GQL']
dealias = 3/2
nz=int(args['--nz'])
ntheta=int(args['--nt'])
nr=int(args['--nr'])
willis=args['--willis']
single_mode = args['--single_mode']
theta_symmetry = args['--theta_symmetry']
perturb_restart = args['--perturb_restart']

if GQL!=None:
    GQL=True
    Lambda_z = int(args['--Lambda_z'])
    Lambda_theta = int(args['--Lambda_theta'])
    logger.info(f"Running with GQL, (Lambda_z, Lambda_theta) = ({Lambda_z}, {Lambda_theta})")
mesh_1 = int(args['--mesh_1'])
mesh_2 = int(args['--mesh_2'])
m1 = int(args['--m'])
restart = args['--restart']
if restart == 'None':
    restart = None
run_note = args['--run_note']
if run_note == 'None':
    run_note = None
stop_time = args['--stop_time']
if stop_time == 'None':
    stop_time = None
else:
    stop_time = float(stop_time)
Ltheta = 2*np.pi
if theta_symmetry:
    logger.info("Running with symmetry restricted theta domain.")
    try:
        Ltheta /= m1
    except ZeroDivisionError:
        raise ZeroDivisionError("m1 is zero. Symmmetry restriction not possible")
if GQL:
    root_folder = "TC_3d_re_{:e}_eta_{:e}_mu_{:e}_Gamma_{:e}_M1_{:d}_{:d}_{:d}_{:d}_GQL_Lambdaz_{:d}_Lambdat_{:d}/".format(Re1,eta,mu,Gamma,m1,nz,ntheta,nr,Lambda_z, Lambda_theta)
else:
    root_folder = "TC_3d_re_{:e}_eta_{:e}_mu_{:e}_Gamma_{:e}_M1_{:d}_{:d}_{:d}_{:d}/".format(Re1,eta,mu,Gamma,m1,nz,ntheta,nr)
if willis:
    root_folder = root_folder.strip("/")
    root_folder += "_willis/"
elif single_mode:
    root_folder = root_folder.strip("/")
    root_folder += "_single_mode/"
if run_note:
    root_folder = root_folder.strip("/")
    root_folder += "_{}/".format(run_note)
if theta_symmetry:
    root_folder = root_folder.strip("/")
    root_folder += "_theta_symmetry/"
path = 'results/'+root_folder
if rank==0:
    if not os.path.exists(path):
        os.mkdir(path)
#derived parameters
R1 = eta/(1. - eta)
R2 = 1./(1-eta)
Omega1 = 1/R1
Omega2 = mu*Omega1
nu = R1 * Omega1/Re1
midpoint = R1 + (R2-R1) / 2
Lz = Gamma 
#Taylor Number
Ta = ((1+eta)**4/(64*eta**2) ) * ( (R2-R1)**2 * (R2+R1)**2 * (Omega1-Omega2)**2 ) / nu**2 #Grossman lohse
Ro_inv = (2 * Omega2 * (R2-R1) ) /  (np.abs(Omega1-Omega2)*R1 )
logger.info("Re:{:.3e}, eta:{:.4e}, mu:{:.4e}".format(Re1,eta,mu))
logger.info("Taylor Number:{:.2e}, Ro^(-1):{:.2e}".format(Ta,Ro_inv))
logger.info("Lz set to {:.6e}".format(Lz))

variables = ['u','ur','v','vr','w','wr','p']
# domain
z_basis = de.Fourier('z', nz, interval=[0., Lz], dealias=dealias)
theta_basis = de.Fourier('theta', ntheta, interval=[0., Ltheta], dealias=dealias)
r_basis = de.Chebyshev('r', nr, interval=[R1, R2], dealias=dealias)
bases = [z_basis, theta_basis, r_basis]
domain = de.Domain(bases, grid_dtype=np.float64, mesh=[mesh_1, mesh_2])  
# problem
problem = de.IVP(domain, variables=variables)
de.operators.parseables["Project"] = GQLProject
# parameters
problem.parameters['eta'] = eta
problem.parameters['mu'] = mu
problem.parameters['Lz'] = Lz
problem.parameters['nu'] = nu
problem.parameters['R1'] = R1
problem.parameters['R2'] = R2
problem.parameters['pi'] = np.pi
#Substitutions
"""
ASSUMES THE EQUATION SET IS PREMULTIPLIED BY A POWER OF r (SEE BELOW)!!!
Lap_s --> scalar laplacian
Lap_r --> r component of vector laplacian
Lap_t --> theta component of vector laplacian
Lap_z --> z component of vector laplacian
"""
problem.substitutions['A'] = '(1/eta - 1.)*(mu-eta**2)/(1-eta**2)'
problem.substitutions['B'] = 'eta*(1-mu)/((1-eta)*(1-eta**2))'
problem.substitutions['v0'] = 'A*r + B/r'       #background profile
problem.substitutions['dv0dr'] = 'A - B/(r*r)'  #d/dr of background forcing
problem.substitutions['v_tot'] = 'v0 + v'       #total velocity in v direction. (azimuthal)
problem.substitutions['vel_sum_sq'] = 'u**2 + v_tot**2 + w**2'
problem.substitutions['plane_avg_r(A)'] = 'integ(integ(r*A, "z"),"theta")/(2*pi*r*Lz)'
#problem.substitutions['plane_avg_z(A)'] = 'integ(integ(A, "r"),"theta")/Lz'
problem.substitutions['vol_avg(A)']   = 'integ(r*A)/(pi*(R2**2 - R1**2)*Lz)'
problem.substitutions['probe(A)'] = 'interp(A,r={}, theta={}, z={})'.format(R1 + 0.5, 0., Lz/2.)
problem.substitutions['KE'] = '0.5*vel_sum_sq'
problem.substitutions['perturb_KE'] = '0.5*(u**2 + v**2 + w**2)'
problem.substitutions['u_rms'] = 'sqrt(u*u)'
problem.substitutions['v_rms'] = 'sqrt(v*v)'
problem.substitutions['w_rms'] = 'sqrt(w*w)'
problem.substitutions['Re_rms'] = 'sqrt(vel_sum_sq)*Lz/nu'
problem.substitutions['epicyclic_freq_sq']  = 'dr(r*r*v*v)/(r*r*r)'
problem.substitutions['enstrophy'] = '0.5*((dtheta(w)/r - dz(v_tot))**2 + (dz(u) - wr )**2 + (vr + dv0dr + v_tot/r - dtheta(u))**2)'

 # not pre-multiplied...don't use this in an equation!
problem.substitutions['DivU'] = "ur + u/r + dtheta(v)/r + dz(w)" 

# assume pre-multiplication by r*r
problem.substitutions['Lap_s(f, f_r)'] = "r*r*dr(f_r) + r*f_r + dtheta(dtheta(f)) + r*r*dz(dz(f))"
problem.substitutions['Lap_r'] = "Lap_s(u, ur) - u - 2*dtheta(v)"
problem.substitutions['Lap_t'] = "Lap_s(v, vr) - v + 2*dtheta(u)"
problem.substitutions['Lap_z'] = "Lap_s(w, wr)"
if GQL:
    # substitutions for projecting onto the low and high modes
    problem.substitutions['Project_high(A)'] = "Project(A,[{:d},{:d}],'h')".format(Lambda_z, Lambda_theta)
    problem.substitutions['Project_low(A)'] = "Project(A,[{:d},{:d}],'l')".format(Lambda_z, Lambda_theta)
    # projected velocities
    problem.substitutions['u_l'] = "Project_low(u)"
    problem.substitutions['u_h'] = "Project_high(u)"
    problem.substitutions['v_l'] = "Project_low(v)"
    problem.substitutions['v_h'] = "Project_high(v)"
    problem.substitutions['w_l'] = "Project_low(w)"
    problem.substitutions['w_h'] = "Project_high(w)"
    # r momentum (GQL)
    problem.add_equation("r*r*dt(u) - nu*Lap_r - 2*r*v0*v + r*v0*dtheta(u) + r*r*dr(p) = - Project_high(r*r*u_h*dr(u_l)+ r*v_h*dtheta(u_l) + r*r*w_h*dz(u_l) - r*v_h*v_l + r*r*u_l*dr(u_h) +r*v_l*dtheta(u_h) + r*r*w_l*dz(u_h) - r*v_h*v_l) - Project_low(r*r*u_h*dr(u_h)+ r*v_h*dtheta(u_h) + r*r*w_h*dz(u_h) - r*v_h*v_h + r*r*u_l*dr(u_l) + r*v_l*dtheta(u_l) + r*r*w_l*dz(u_l) - r*v_l*v_l)")
    # theta momentum (GQL)
    problem.add_equation("r*r*dt(v) - nu*Lap_t + r*r*dv0dr*u + r*v0*u + r*v0*dtheta(v) + r*dtheta(p) = - Project_high(r*r*u_h*dr(v_l) + r*v_h*dtheta(v_l) + r*r*w_h*dz(v_l) + r*v_h*u_l + r*r*u_l*dr(v_h) + r*v_l*dtheta(v_h) + r*r*w_l*dz(v_h) + r*v_l*u_h) - Project_low(r*r*u_h*dr(v_h) + r*v_h*dtheta(v_h) + r*r*w_h*dz(v_h) + r*v_h*u_h + r*r*u_l*dr(v_l) + r*v_l*dtheta(v_l) + r*r*w_l*dz(v_l) + r*v_l*u_l)")
    # z momentum (GQL)
    problem.add_equation("r*r*dt(w) - nu*Lap_z + r*r*dz(p) + r*v0*dtheta(w) = - Project_high(r*r*u_h*dr(w_l) + r*v_h*dtheta(w_l) + r*r*w_h*dz(w_l) + r*r*u_l*dr(w_h) + r*v_l*dtheta(w_h) + r*r*w_l*dz(w_h)) - Project_low(r*r*u_h*dr(w_h) + r*v_h*dtheta(w_h) + r*r*w_h*dz(w_h) + r*r*u_l*dr(w_l) + r*v_l*dtheta(w_l) + r*r*w_l*dz(w_l))")
else:
    # DNS nonlinear terms
    problem.substitutions['UdotGrad_s(f, f_r)'] = "r*r*u*f_r + r*v*dtheta(f) + r*r*w*dz(f)"
    problem.substitutions['UdotGrad_r'] = "UdotGrad_s(u, ur) - r*v*v"
    problem.substitutions['UdotGrad_t'] = "UdotGrad_s(v, vr) + r*u*v"
    problem.substitutions['UdotGrad_z'] = "UdotGrad_s(w, wr)"
    # momentum equations
    problem.add_equation("r*r*dt(u) - nu*Lap_r - 2*r*v0*v + r*v0*dtheta(u) + r*r*dr(p) = - UdotGrad_r")
    problem.add_equation("r*r*dt(v) - nu*Lap_t + r*r*dv0dr*u + r*v0*u + r*v0*dtheta(v) + r*dtheta(p)  = -UdotGrad_t  ")
    problem.add_equation("r*r*dt(w) - nu*Lap_z + r*r*dz(p) + r*v0*dtheta(w) = -UdotGrad_z")

#continuity
problem.add_equation("r*ur + u + dtheta(v) + r*dz(w) = 0")
#Auxillilary equations
problem.add_equation("ur - dr(u) = 0")
problem.add_equation("vr - dr(v) = 0")
problem.add_equation("wr - dr(w) = 0")
#Boundary Conditions
problem.add_bc("left(u) = 0")
problem.add_bc("right(u) = 0", condition="ntheta != 0 or nz != 0")
problem.add_bc("right(p) = 0", condition="ntheta == 0 and nz == 0") #??
problem.add_bc("left(v) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")
#Create solver
solver = problem.build_solver(de.timesteppers.RK443)
#intial conditions
u = solver.state['u']
ur = solver.state['ur']
v = solver.state['v']
vr = solver.state['vr']
w = solver.state['w']
wr = solver.state['wr']
r = problem.domain.grid(-1, scales = problem.domain.dealias)
z = problem.domain.grid(0, scales = problem.domain.dealias)
theta = problem.domain.grid(1, scales = problem.domain.dealias)
A0 = 1e-3
write_mode = "overwrite"
if restart:
    write_mode = "append"
    logger.info("Restarting from file {}".format(restart))
    write, last_dt = solver.load_state(restart, -1)
    if perturb_restart:
        logger.info("perturbing restarted velocity field with A0 = {}".format(A0))
        Ar, Atheta, Az = random_vector_potential(domain, R1, R2)
        for vel in [u, v, w]:
            vel.set_scales(domain.dealias, keep_data=True)
        # Curl of A
        u['g'] += A0 * (Az.differentiate('theta')['g']/r - Atheta.differentiate('z')['g'])
        v['g'] += A0 * (Ar.differentiate('z')['g'] - Az.differentiate('r')['g'])
        w['g'] += A0 * (Atheta['g']/r + Atheta.differentiate('r')['g'] - Ar.differentiate('theta')['g']/r)
        u.differentiate('r', out=ur)
        v.differentiate('r', out=vr)
        w.differentiate('r', out=wr)
elif willis:
    ## Willis & Bahrenghi ICs
    logger.info("Using initial conditions from Willis's PhD thesis")
    u.set_scales(domain.dealias, keep_data=False)
    w.set_scales(domain.dealias, keep_data=False)
    x = r - R1
    kz = 2*np.pi/Lz
    logger.info('kz : {}'.format(kz))
    u['g'] = A0 * kz**2 * x**2 * (1-x)**2 * np.sin(kz*z)
    w['g'] = A0 * (kz * x**2 * (1-x)**2 * np.cos(kz*z)/r + 2*kz*np.cos(kz*z) * ((1-x)**2 * x - x**2 * (1 - x)) - (x**2 * (1 - x)**2)/r * m1 * np.cos(m1*theta))
    u.differentiate('r',out=ur)
    w.differentiate('r',out=wr)
elif single_mode:
    logger.info("using single-mode initial conditions with kz = 1, m = {}; amplitude {}".format(m1, A0))
    kz = 2*np.pi/Lz
    for vel in [v, w]:
        vel.set_scales(domain.dealias, keep_data=False)

    v['g'] = A0 * np.sin(np.pi*(r-R1)) * np.sin(kz*z)
    w['g'] = A0 * np.sin(np.pi*(r-R1)) * np.sin(m1*theta)
    v.differentiate('r', out=vr)
    w.differentiate('r', out=wr)
else:
    logger.info("Using incompressible noise initial conditions in (u, v, w) with amplitude A0 = {}.".format(A0))
    Ar, Atheta, Az = random_vector_potential(domain, R1, R2)
    for vel in [u, v, w]:
        vel.set_scales(domain.dealias, keep_data=True)
    # Curl of A
    u['g'] = A0 * (Az.differentiate('theta')['g']/r - Atheta.differentiate('z')['g'])
    v['g'] = A0 * (Ar.differentiate('z')['g'] - Az.differentiate('r')['g'])
    w['g'] = A0 * (Atheta['g']/r + Atheta.differentiate('r')['g'] - Ar.differentiate('theta')['g']/r)
    u.differentiate('r', out=ur)
    v.differentiate('r', out=vr)
    w.differentiate('r', out=wr)
#Setting Simulation Runtime
omega1 = 1/eta - 1.
period = 2*np.pi/omega1
t_visc = Re1
if stop_time:
    logger.info(f"running for {stop_time} time units")
    solver.stop_sim_time = solver.sim_time + stop_time
else:
    logger.info(f"running for one viscous time: t_visc = {t_visc}")
    solver.stop_sim_time = solver.sim_time + t_visc
solver.stop_wall_time = 428400. # 5 days - 1 hour
solver.stop_iteration = np.inf #2000
#CFL stuff
CFL = flow_tools.CFL(solver, initial_dt=1e-2, cadence=5, safety=0.3, max_change=1.5, min_change=0.5, max_dt=0.1, threshold=0.1)
CFL.add_velocities(('u', 'v', 'w'))
# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("abs(DivU)", name='divu')
flow.add_property("integ(r*KE)", name='KE')
flow.add_property("integ(r*enstrophy)", name='enstrophy')
dt = CFL.compute_dt()

geo_factor = 1

#Analysis
output_time_cadence = 0.1*period
scalar_output_time_cadence = output_time_cadence/10.
analysis_tasks = []
logger.info("path= {}".format(path))
# checkpoints
snapshots = solver.evaluator.add_file_handler(path + 'snapshots',sim_dt=period,max_writes=1, mode=write_mode)
snapshots.add_system(solver.state)
analysis_tasks.append(snapshots)
# slices
analysis_slice = solver.evaluator.add_file_handler(path+"/slices", parallel=False, sim_dt=output_time_cadence, mode=write_mode)
analysis_slice.add_task("interp(u,r={})".format(midpoint), name="u_r_slice",scales=4)
analysis_slice.add_task("interp(v,r={})".format(midpoint), name="v_r_slice",scales=4)
analysis_slice.add_task("interp(w,r={})".format(midpoint), name="w_r_slice",scales=4)
analysis_slice.add_task("interp(u,theta=0)", name="u_t_slice",scales=4)
analysis_slice.add_task("interp(v,theta=0)", name="v_t_slice",scales=4)
analysis_slice.add_task("interp(w,theta=0)", name="w_t_slice",scales=4)
analysis_slice.add_task("interp(u,z={})".format(Lz/2), name="u_z_slice",scales=4)
analysis_slice.add_task("interp(v,z={})".format(Lz/2), name="v_z_slice",scales=4)
analysis_slice.add_task("interp(w,z={})".format(Lz/2), name="w_z_slice",scales=4)
analysis_tasks.append(analysis_slice)
# profiles
analysis_profile = solver.evaluator.add_file_handler(path+"/profiles", max_writes=20, parallel=False, mode=write_mode)
analysis_slice.add_task("plane_avg_r(v_tot)", name="v_tot")
analysis_slice.add_task("plane_avg_r(u_rms)", name="u_rms")
analysis_slice.add_task("plane_avg_r(v_rms)", name="v_rms")
analysis_slice.add_task("plane_avg_r(w_rms)", name="w_rms")
analysis_slice.add_task("plane_avg_r(Re_rms)", name="Re_rms")
analysis_slice.add_task("plane_avg_r(epicyclic_freq_sq)", name="epicyclic_freq_sq")
analysis_tasks.append(analysis_profile)
# spectra
analysis_spectra = solver.evaluator.add_file_handler(path+"/spectra", max_writes=20, parallel=False, sim_dt=output_time_cadence, mode=write_mode)
analysis_spectra.add_task("interp(u, r={})".format(midpoint), name="uc", layout="c")
analysis_spectra.add_task("interp(v, r={})".format(midpoint), name="vc", layout="c")
analysis_spectra.add_task("interp(w, r={})".format(midpoint), name="wc", layout="c")
analysis_tasks.append(analysis_spectra)
# scalars
analysis_scalar = solver.evaluator.add_file_handler(path+"/scalar", parallel=False,sim_dt=scalar_output_time_cadence, mode=write_mode)
analysis_scalar.add_task("integ(r*KE)", name="KE")
analysis_scalar.add_task("vol_avg(u_rms)", name="u_rms")
analysis_scalar.add_task("vol_avg(v_rms)", name="v_rms")
analysis_scalar.add_task("vol_avg(w_rms)", name="w_rms")
analysis_scalar.add_task("vol_avg(Re_rms)", name="Re_rms")
analysis_scalar.add_task("probe(u)", name="u_probe")
analysis_scalar.add_task("probe(w)", name="w_probe")
analysis_scalar.add_task("integ(r*enstrophy)", name="enstrophy")
analysis_scalar.add_task("integ(r*perturb_KE)", name="pertubation_KE")
analysis_tasks.append(analysis_scalar)

logger.info("Starting main loop...")
cond_number(solver)
while solver.ok:
    solver.step(dt)
    if (solver.iteration-1) % 100 == 0:
        logger.info('Iteration: %i, Time: %e, t/P_in: %e, t/t_visc: %e, dt: %e' %(solver.iteration, solver.sim_time, solver.sim_time/period, solver.sim_time/Re1, dt))
        logger.info('Max |divu| = {}'.format(flow.max('divu')))
        logger.info('Total KE per Lz = {}'.format(geo_factor*flow.max('KE')/Lz))
        logger.info('Total enstrophy per Lz = {}'.format(geo_factor*flow.max('enstrophy')/Lz))
        # Hermitian projection
        for field in solver.state.fields:
            field.require_grid_space()

    dt = CFL.compute_dt()
solver.evaluate_handlers_now(dt)
solver.log_stats()
for task in analysis_tasks:
    post.merge_analysis(task.base_path)

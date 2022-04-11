"""
Usage:
  taylor_couette_3d.py --re=<reynolds> --eta=<eta> --m=<initial_m> [--mu=<mu>] [--ar=<Gamma>] [--restart=<restart>] [--restart_file=<restart_file>] --mesh_1=<mesh_1> --mesh_2=<mesh_2> [--GQL=<GQL>] [--Lambda_z=<Lambda_z>] [--Lambda_theta=<Lambda_theta>] [--willis] [--single_mode] [--run_note=<run_note>] [--theta_symmetry]
  taylor_couette_3d.py

Options:
  --re=<reynolds>  Reynolds number for simulation
  --eta=<eta>      Eta - ratio of R1/R2
  --m=<initial_m>  M1 mode to begin initial conditions
  --mu=<mu>        mu = Omega2/Omega1 [default: 0]
  --ar=<Gamma>     Aspect ratio (height/width)
  --mesh_1=<mesh_1> First mesh core count
  --mesh_2=<mesh_2> Second mesh core count
  --restart=<restart> True or False
  --restart_file=<restart_file> Point to a merged snapshots_s1.h5 file
  --GQL=<GQL> True or False
  --Lambda_z=<Lambda_z> Specify an integer cutoff to seperate low and high modes for z
  --Lambda_theta=<Lambda_theta> Specify an integer cutoff to seperate low and high modes for theta
  --willis  Use Willis ICs [default: False]  
  --single_mode  Use single mode ICs [default: False]
  --run_note=<run_note>  Note to add to run directory name [default: None]
  --theta_symmetry  Restrict theta to 2pi/m1 [default: False]
"""

import numpy as np
import h5py
from dedalus.extras import flow_tools
from dedalus import public as de
import time
import logging
from docopt import docopt
import os
import subprocess
from mpi4py import MPI
from GQLProjection import GQLProjection
from filter_field import filter_field

def cond_number(solver):
    cond_nums = np.zeros(len(solver.pencils))
    for i,p in enumerate(solver.pencils):
        cond_nums[i] = np.linalg.cond(p.L.todense())
    comm = solver.domain.dist.comm
    global_max = np.array([0.,])
    local_max = np.array([cond_nums.max(),])
    comm.Reduce(local_max,global_max, MPI.MAX)
    if comm.rank == 0:
        logger.info("Maximum condition number is {:e}".format(global_max[0]))


logger = logging.getLogger(__name__)

de.operators.parseables["Project"] = GQLProjection

comm=MPI.COMM_WORLD
rank=comm.Get_rank()

args=docopt(__doc__)
Re1=float(args['--re'])
eta=np.float(args['--eta'])
mu = np.float(args['--mu'])
Gamma = float(args['--ar'])
GQL = args['--GQL']

willis=args['--willis']
single_mode = args['--single_mode']
theta_symmetry = args['--theta_symmetry']

if GQL!=None:
    GQL=True
    Lambda_z = int(args['--Lambda_z'])
    Lambda_theta = int(args['--Lambda_theta'])
mesh_1 = int(args['--mesh_1'])
mesh_2 = int(args['--mesh_2'])
m1 = int(args['--m'])
restart = bool(args['--restart'])
run_note = args['--run_note']
if run_note == 'None':
    run_note = None

Ltheta = 2*np.pi
if theta_symmetry:
    logger.info("Running with symmetry restricted theta domain.")
    try:
        Ltheta /= m1
    except ZeroDivisionError:
        raise ZeroDivisionError("m1 is zero. Symmmetry restriction not possible")

if restart==True:
    restart_file = str(args['--restart_file'])
"""
delta = R2 - R1
mu = Omega2/Omega1
eta = R1/R2

scale [L] = delta
scale [T] = delta/(R1 Omega1)
scale [V] = R1 Omega1

Default parameters from Barenghi (1991, J. Comp. Phys.).
"""

#eta = 0.8770
# ~ Re1 = 80
#Lz = 2.0074074463832545
Sc = 1
dealias = 3/2
nz=64
ntheta=128
nr=32

#eta_string = "{:.4e}".format(eta).replace(".","-")


if GQL:
    root_folder = "TC_3d_re_{:e}_eta_{:e}_Gamma_{:e}_M1_{:d}_{:d}_{:d}_{:d}_GQL_Lambdaz_{:d}_Lambdat_{:d}/".format(Re1,eta,Gamma,m1,nz,ntheta,nr,Lambda_z, Lambda_theta)
else:
    root_folder = "TC_3d_re_{:e}_eta_{:e}_Gamma_{:e}_M1_{:d}_{:d}_{:d}_{:d}/".format(Re1,eta,Gamma,m1,nz,ntheta,nr)

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
        # os.chdir(path) # Allows us to run slurm_analysis in run_GQL.sh or run_DNS.sh
    elif restart==False:
        logger.info('Folder for run already exists.')
        logger.info('Use restart, rename existing folder, or change parameters')
        #subprocess.call(['analysis_scripts/./kill_script.sh'])
sim_name=path
restart=False

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
Ta1 = (2*Omega1**2*(R2-R1)**4*eta**2 ) / (nu**2 *(11-eta**2)  )
Ro_inv = (2 * Omega2 * (R2-R1) ) /  (np.abs(Omega1-Omega2)*R1 )

logger.info("Re:{:.3e}, eta:{:.4e}, mu:{:.4e}".format(Re1,eta,mu))
logger.info("Taylor Number:{:.2e}, Ro^(-1):{:.2e}".format(Ta,Ro_inv))

logger.info("Lz set to {:.6e}".format(Lz))

variables = ['u','ur','v','vr','w','wr','p']

#domain
z_basis = de.Fourier('z', nz, interval=[0., Lz], dealias=dealias)
theta_basis = de.Fourier('theta', ntheta, interval=[0., Ltheta], dealias=dealias)
r_basis = de.Chebyshev('r', nr, interval=[R1, R2], dealias=dealias)

bases = [z_basis, theta_basis, r_basis]
# ~ bases = t_bases + r_basis
domain = de.Domain(bases, grid_dtype=np.float64, mesh=[mesh_1,mesh_2])  

#problem
problem = de.IVP(domain, variables=variables)

#params into equations
problem.parameters['eta'] = eta
problem.parameters['mu'] = mu
problem.parameters['Lz'] = Lz
problem.parameters['nu'] = nu
problem.parameters['R1'] = R1
problem.parameters['R2'] = R2
problem.parameters['pi'] = np.pi

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

# substitutions
if GQL:

    # substitutions for projecting onto the low and high modes
    problem.substitutions['Project_high(A)'] = "Project(A,[{:d},{:d}],'h')".format(Lambda_z, Lambda_theta)
    problem.substitutions['Project_low(A)'] = "Project(A,[{:d},{:d}],'l')".format(Lambda_z, Lambda_theta)

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
r_in = R1


A0 = 1e-3
if restart == True:
	logger.info("Restarting from file {}".format(restart_file))
	write, last_dt = solver.load_state(restart_file, -1)
elif willis:
    ## Willis & Bahrenghi ICs
    logger.info("Using initial conditions from Willis's PhD thesis")


    #m1=3
    u.set_scales(domain.dealias, keep_data=False)
    w.set_scales(domain.dealias, keep_data=False)
    x = r - r_in
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

    v['g'] = A0 * np.sin(np.pi*(r-r_in)) * np.sin(kz*z)
    w['g'] = A0 * np.sin(np.pi*(r-r_in)) * np.sin(m1*theta)
    v.differentiate('r', out=vr)
    w.differentiate('r', out=wr)
else:
    # Random perturbations to v in (r, z)
    gshape = domain.dist.grid_layout.global_shape(scales=domain.dealias)
    slices = domain.dist.grid_layout.slices(scales=domain.dealias)
    rand = np.random.RandomState(seed=42)

    logger.info("Using incompressible noise initial conditions in (u, v, w) with amplitude A0 = {}.".format(A0))
    filter_fraction = 0.5
    Ar = domain.new_field()
    Atheta = domain.new_field()
    Az = domain.new_field()

    fr= (4*(r-R1)*(R2-r))**4
    for A in [Ar, Atheta, Az]:
        A.set_scales(domain.dealias, keep_data=False)
        A['g'] = rand.standard_normal(gshape)[slices]
        A.set_scales(0.5,keep_data=True)
        A['c']
        A['g']
        A.set_scales(domain.dealias,keep_data=True)
        A['g'] *= fr
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
solver.stop_sim_time = 10 * period # np.inf #6*period
solver.stop_wall_time = 24 * 3600. #np.inf # This is in seconds
solver.stop_iteration = np.inf #2000

#CFL stuff
CFL = flow_tools.CFL(solver, initial_dt=1e-2, cadence=5, safety=0.3, max_change=1.5, min_change=0.5, max_dt=1, threshold=0.1)
CFL.add_velocities(('u', 'v', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("abs(DivU)", name='divu')
flow.add_property("integ(r*KE)", name='KE')
flow.add_property("integ(r*enstrophy)", name='enstrophy')

dt = CFL.compute_dt()
# Main loop


geo_factor = 1

#Analysis

output_time_cadence = 0.1*period
scalar_output_time_cadence = output_time_cadence/100.

# ~ analysis = solver.evaluator.add_file_handler('taylor_couette',scalar_output_time_cadence,max_writes=np.inf)
logger.info("sim_name= {}".format(sim_name))
snapshots = solver.evaluator.add_file_handler(sim_name + 'snapshots',sim_dt=output_time_cadence,max_writes=np.inf)
snapshots.add_system(solver.state)

#Analysis files
Jeffs_analysis=True
if Jeffs_analysis:
    analysis_slice = solver.evaluator.add_file_handler(sim_name+"/slices", parallel=False, sim_dt=output_time_cadence)
    analysis_slice.add_task("interp(u,r={})".format(midpoint), name="u_slice",scales=4)
    analysis_slice.add_task("interp(v,r={})".format(midpoint), name="v_slice",scales=4)
    analysis_slice.add_task("interp(w,r={})".format(midpoint), name="w_slice",scales=4)

    analysis_slice.add_task("interp(KE, z=0)", name="KE")
    analysis_slice.add_task("plane_avg_r(v_tot)", name="v_tot")
    analysis_slice.add_task("plane_avg_r(u_rms)", name="u_rms")
    analysis_slice.add_task("plane_avg_r(v_rms)", name="v_rms")
    analysis_slice.add_task("plane_avg_r(w_rms)", name="w_rms")
    analysis_slice.add_task("plane_avg_r(Re_rms)", name="Re_rms")
    analysis_slice.add_task("plane_avg_r(epicyclic_freq_sq)", name="epicyclic_freq_sq")
    analysis_slice.add_task("integ(r*v, 'z')", name='Angular Momentum')
    
    analysis_profile = solver.evaluator.add_file_handler(sim_name+"/profiles", max_writes=20, parallel=False)

    analysis_spectra = solver.evaluator.add_file_handler(sim_name+"/spectra", max_writes=20, parallel=False, sim_dt=output_time_cadence)
    analysis_spectra.add_task("interp(u, r={})".format(midpoint), name="uc", layout="c")
    analysis_spectra.add_task("interp(v, r={})".format(midpoint), name="vc", layout="c")
    analysis_spectra.add_task("interp(w, r={})".format(midpoint), name="wc", layout="c")
    
    analysis_scalar = solver.evaluator.add_file_handler(sim_name+"/scalar", parallel=False,sim_dt=scalar_output_time_cadence)
    analysis_scalar.add_task("integ(r*KE)", name="KE")
    analysis_scalar.add_task("vol_avg(u_rms)", name="u_rms")
    analysis_scalar.add_task("vol_avg(v_rms)", name="v_rms")
    analysis_scalar.add_task("vol_avg(w_rms)", name="w_rms")
    analysis_scalar.add_task("vol_avg(Re_rms)", name="Re_rms")
    analysis_scalar.add_task("probe(w)", name="w_probe")
    analysis_scalar.add_task("integ(r*enstrophy)", name="enstrophy")
    analysis_scalar.add_task("integ(r*perturb_KE)", name="pertubation_KE")

logger.info("Starting main loop...")
cond_number(solver)
start_time = time.time()
while solver.ok:
    solver.step(dt)
    if (solver.iteration-1) % 100 == 0:
        logger.info('Iteration: %i, Time: %e, Inner rotation periods: %e, dt: %e' %(solver.iteration, solver.sim_time, solver.sim_time/period, dt))
        logger.info('Max |divu| = {}'.format(flow.max('divu')))
        logger.info('Total KE per Lz = {}'.format(geo_factor*flow.max('KE')/Lz))
        logger.info('Total enstrophy per Lz = {}'.format(geo_factor*flow.max('enstrophy')/Lz))

        # Hermitian projection
        for field in solver.state.fields:
            field.require_grid_space()

    dt = CFL.compute_dt()

end_time = time.time()

logger.info('Simulation run time: {:f}'.format(end_time-start_time))

logger.info('Time per iteration: {:f}'.format((end_time-start_time)/solver.iteration))

logger.info('Simulation ended')









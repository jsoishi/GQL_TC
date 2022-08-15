from pathlib import Path
import h5py
from series import RunSeries
import matplotlib.pyplot as plt
import numpy as np
import paths
plt.style.use('prl')
#path = RunSeries(paths.rei640_reo_1359, 'rei640_reo_-1359')
path = RunSeries(paths.rei900_reo_3398, 'rei900_reo_-3398')

t_cutoff = 1.7
filetype = Path("slices/slices.h5")
path.data = []
path.r = []
path.sim_time = []
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
for i, r in enumerate(path.run_list[:-2]):
    filename = path.root_dir/Path(r)/filetype
    with h5py.File(filename, "r") as df:
        path.data.append(df['tasks/v_tot'][:,0,0,:])
        path.r.append(df['scales/r/1.0'][:])
        path.sim_time.append(df['scales/sim_time'][:])
    label = path.Lambdaz[i]
    eta = path.eta[i]
    mu = path.mu[i]
    t_visc = 1/path.Re[i]
    A = (1/eta - 1.)*(mu-eta**2)/(1-eta**2)
    B = eta*(1-mu)/((1-eta)*(1-eta**2))
    r = path.r[i]
    vr = A*r + B/r
    dns_profile = path.data[0][path.sim_time[0]/t_visc>t_cutoff,:].mean(axis=0)
    r_profile = path.data[i][path.sim_time[i]/t_visc>t_cutoff,:].mean(axis=0)
    if label is None:
        label = 'DNS'
        c = 'k'
        lw = 3
    else:
        label = f"$\Lambda = {label}$"
        c = colors[i]
        lw = 2
        l2_norm = ((r_profile - dns_profile)**2).mean()
        print(f"Lambda = {label} l2_norm = {l2_norm}")
    print(f"Lambda = {label} number of samples = {path.data[i].shape[0]}")
    plt.plot(r, r_profile -vr, label=label, color=c, lw=lw)
plt.xlabel(r"$r/d$")
plt.ylabel(r"$\left< v' \right>$")
plt.xlim(path.r[0].min(),path.r[0].max())
plt.legend(loc=(0.5,0.05),fontsize=16)
plt.tight_layout()
plt.savefig(f"../figs/{path.name}_vmean_profile.pdf")

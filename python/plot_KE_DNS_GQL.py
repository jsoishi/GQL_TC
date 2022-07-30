from pathlib import Path
import h5py
from series import RunSeries
import matplotlib.pyplot as plt
import numpy as np
#from paths import rei640_reo_1359
from paths import rei900_reo_3398

plt.style.use('prl')
#path = RunSeries(rei640_reo_1359, 'rei640_reo_-1359')
path = RunSeries(rei900_reo_3398, 'rei900_reo_-3398')

filetype = Path("scalar/scalar.h5")
path.data = []
path.t = []
t_visc = path.Re[0]
print(f"t_visc = {t_visc}")
for i, r in enumerate(path.run_list):
    filename = path.root_dir/Path(r)/filetype
    with h5py.File(filename, "r") as df:
        path.data.append(df['tasks/pertubation_KE'][:,0,0,0])
        path.t.append(df['scales/sim_time'][:])
    label = path.Lambdaz[i]

    if label is None:
        label = '\infty'
    print(f"Lambda = {label} number of samples = {path.data[i].shape[0]}")
    plt.plot(path.t[i]/t_visc, path.data[i], label=f"$\Lambda = {label}$")
plt.xlabel(r"$t/\tau_\nu$")
plt.ylabel(r"$E_{kin}$")
plt.xlim(1,2.5)
plt.legend()
plt.tight_layout()
plt.savefig(f"../figs/{path.name}_KE_vs_t.pdf")

from pathlib import Path
import h5py
from series import RunSeries
import matplotlib.pyplot as plt
import numpy as np
from paths import rei640_reo_1359

plt.style.use('prl')
path = RunSeries(rei640_reo_1359, 'rei640_reo_-1359')

filetype = Path("slices/slices_s1.h5")
path.data = []
path.r = []

for i, r in enumerate(path.run_list):
    filename = path.root_dir/Path(r)/filetype
    with h5py.File(filename, "r") as df:
        path.data.append(df['tasks/v_tot'][:,0,0,:])
        path.r.append(df['scales/r/1.0'][:])
    label = path.Lambdaz[i]

    if label is None:
        label = '\infty'
    print(f"Lambda = {label} number of samples = {path.data[i].shape[0]}")
    plt.plot(path.r[i],path.data[i][-25:,:].mean(axis=0), label=f"$\Lambda = {label}$")
plt.xlabel(r"$r/d$")
plt.ylabel(r"$\left< v \right>$")
plt.xlim(path.r[0].min(),path.r[0].max())
plt.legend()
plt.tight_layout()
plt.savefig(f"../figs/{path.name}_vmean_profile.pdf")

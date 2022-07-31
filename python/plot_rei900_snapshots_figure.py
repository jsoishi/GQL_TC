from pathlib import Path
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import paths
from series import RunSeries
from dedalus.extras import plot_tools
plt.style.use('prl')

#path = RunSeries(paths.rei640_reo_1359, 'rei640_reo_-1359')
path = RunSeries(paths.rei900_reo_3398, 'rei900_reo_-3398')

aspect=0.9 #np.sqrt(2)
w = 5
h = aspect*w

def load_data(df, Rmid, index=-1):
    dset = np.sqrt(df['tasks']['u_r_slice'][index, :,:,0]**2 + df['tasks']['v_r_slice'][index, :,:,0]**2 + df['tasks']['w_r_slice'][index, :,:,0]**2)

    image_axes = [2, 1]
    data_slices = [index, slice(None), slice(None),0]
    xmesh, ymesh, vdata = plot_tools.get_plane(df['tasks/u_r_slice'], 2, 1, data_slices)
    xmesh *= Rmid

    return xmesh, ymesh, dset


fig = plt.figure(figsize=(w, h))
grid = AxesGrid(fig, 111,
                nrows_ncols=(3, 2),
                axes_pad=0.05,
                cbar_mode='single',
                cbar_location='top',
                cbar_pad=0.
                )


for i,ax in enumerate(grid):
    if i != 4:
        ax.set_axis_off()
    else:
        ax.set_xlabel(r"$R_{mid}\theta$", fontsize=14)
        ax.set_ylabel("z", fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)
    r = path.run_list[i]
    eta = path.eta[i]
    print(f"plotting {r}")

    filetype = 'slices/slices_s2.h5'
    if not (path.root_dir/Path(r)/filetype).exists():
        filetype = 'slices/slices_s1.h5'
    filename = path.root_dir/Path(r)/filetype

    Rmid = eta/(1-eta) + 0.5
    with h5py.File(filename, "r") as df:
        xgrid, ygrid, data = load_data(df, Rmid)

    im = ax.pcolormesh(xgrid, ygrid, data, vmin=0, vmax=2, rasterized=True)
    label = path.Lambdaz[i]
    if label is None:
        label = '\infty'
    ax.text(2,25,f"$\Lambda = {label}$", color='white')

# when cbar_mode is 'single', for ax in grid, ax.cax = grid.cbar_axes[0]

cbar = ax.cax.colorbar(im)
cbar = grid.cbar_axes[0].colorbar(im)
ax.cax.xaxis.set_label_position('top')
ax.cax.xaxis.set_ticks_position('top')
ax.cax.tick_params(labelsize=12)
cbar.set_label(r"$u'_{rms}$")

plt.savefig("../figs/rei900_snapshots.pdf", dpi=150)

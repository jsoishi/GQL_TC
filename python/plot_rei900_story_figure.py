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


aspect=0.7 #np.sqrt(2)
w = 10
h = aspect*w

def load_data(df, Rmid, index=-1):
    dset = np.sqrt(df['tasks']['u_r_slice'][index, :,:,0]**2 + df['tasks']['v_r_slice'][index, :,:,0]**2 + df['tasks']['w_r_slice'][index, :,:,0]**2)

    image_axes = [2, 1]
    data_slices = [index, slice(None), slice(None),0]
    xmesh, ymesh, vdata = plot_tools.get_plane(df['tasks/u_r_slice'], 2, 1, data_slices)
    xmesh *= Rmid

    return xmesh, ymesh, dset

rect = [0.5,0.15,0.4,0.8]
fig = plt.figure(figsize=(w, h))
ekin_ax = fig.add_axes([0.1,0.15,0.4,0.8])
grid = AxesGrid(fig, rect,
                nrows_ncols=(3, 1),
                axes_pad=0.05,
                # cbar_mode='single',
                # cbar_location='top',
                # cbar_pad=0.
                )
filetype = Path("scalar/scalar.h5")
path.data = []
path.t = []
t_visc = path.Re[0]
print(f"t_visc = {t_visc}")
path_i = 4
indices = [81,89,-4]

r = path.run_list[path_i]
eta = path.eta[path_i]
label = path.Lambdaz[path_i]
filename = path.root_dir/Path(r)/filetype
with h5py.File(filename, "r") as df:
    path.data.append(df['tasks/pertubation_KE'][:,0,0,0])
    path.t.append(df['scales/sim_time'][:])


if label is None:
    label = '\infty'

ekin_ax.plot(path.t[0]/t_visc, path.data[0], label=f"$\Lambda = {label}$")
ekin_ax.set_xlabel(r"$t/\tau_\nu$")
ekin_ax.set_ylabel(r"$E_{kin}$")
ekin_ax.set_xlim(1.8,2.0)
ekin_ax.set_ylim(0,350)
filetype = 'slices/slices_s2.h5'
if not (path.root_dir/Path(r)/filetype).exists():
    filetype = 'slices/slices_s1.h5'
filename = path.root_dir/Path(r)/filetype
print(f"images from {filename}")
with h5py.File(filename, "r") as df:
    for i,ax in enumerate(grid):
        if i != 2:
            ax.set_axis_off()
        else:
            ax.set_xlabel(r"$R_{mid}\theta$", fontsize=14)
            ax.set_ylabel("z", fontsize=14)
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()
            ax.tick_params(axis='both', which='major', labelsize=14)

        print(f"plotting index {indices[i]}")
        Rmid = eta/(1-eta) + 0.5
        xgrid, ygrid, data = load_data(df, Rmid, index=indices[i])

        im = ax.pcolormesh(xgrid, ygrid, data, rasterized=True)
        ekin_ax.axvline(df['scales/sim_time'][indices[i]]/t_visc, alpha=0.4)
        ax.text(2,25,f"{i+1}", color='white', fontsize=14, bbox=dict(boxstyle='circle',facecolor='gray'))
        ekin_ax.text(df['scales/sim_time'][indices[i]]/t_visc-0.0125,320,f"{i+1}", color='white', fontsize=14, bbox=dict(boxstyle='circle',facecolor='gray'))

# when cbar_mode is 'single', for ax in grid, ax.cax = grid.cbar_axes[0]

# cbar = ax.cax.colorbar(im)
# cbar = grid.cbar_axes[0].colorbar(im)
# ax.cax.xaxis.set_label_position('top')
# ax.cax.xaxis.set_ticks_position('top')
# ax.cax.tick_params(labelsize=12)
# cbar.set_label(r"$u'_{rms}$")

plt.savefig("../figs/rei900_story.pdf", dpi=150)

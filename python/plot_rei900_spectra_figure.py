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

index = -1 # time index
eps = 1e-16
aspect=1.6 #np.sqrt(2)
w = 5
h = aspect*w
nk = 100
nm = 100

def load_data(df, Rmid, index=-1):
    dset = df['tasks']['v_r_slice'][index, :,:,0]

    image_axes = [2, 1]
    data_slices = [index, slice(None), slice(None),0]
    xmesh, ymesh, vdata = plot_tools.get_plane(df['tasks/u_r_slice'], 2, 1, data_slices)
    xmesh *= Rmid

    return xmesh, ymesh, dset
def load_snapshot_data(df, Rmid, index=-1):
    mid = int(df['tasks']['v'].shape[-1]/2)
    dset = df['tasks']['v'][index, :,:,mid]

    image_axes = [2, 1]
    data_slices = [index, slice(None), slice(None),mid]
    xmesh, ymesh, vdata = plot_tools.get_plane(df['tasks/v'], 2, 1, data_slices)
    xmesh *= Rmid

    return xmesh, ymesh, dset

def load_data(df, Rmid, index=-1):
    dset = df['tasks']['v_r_slice'][index, :,:,0]

    image_axes = [2, 1]
    data_slices = [index, slice(None), slice(None),0]
    xmesh, ymesh, vdata = plot_tools.get_plane(df['tasks/u_r_slice'], 2, 1, data_slices)
    xmesh *= Rmid

    return xmesh, ymesh, dset


fig = plt.figure(figsize=(w, h))
grid = AxesGrid(fig, 111,
                nrows_ncols=(4, 2),
                axes_pad=0.05,
                cbar_mode='single',
                cbar_location='top',
                cbar_pad=0.
                )


for i,ax in enumerate(grid):
    if i != 6:
        ax.set_axis_off()
    else:
        ax.set_xlabel(r"$m$", fontsize=14)
        ax.set_ylabel(r"$k_z$", fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)
    r = path.run_list[i]
    eta = path.eta[i]
    #print(f"plotting {r}")

    filetype = 'slices/slices_s2.h5'
    if not (path.root_dir/Path(r)/filetype).exists():
        filetype = 'slices/slices_s1.h5'
    if not (path.root_dir/Path(r)/filetype).exists():
        filetype = 'snapshots/snapshots_s4.h5'
    filename = path.root_dir/Path(r)/filetype

    Rmid = eta/(1-eta) + 0.5
    with h5py.File(filename, "r") as df:
        try:
            xgrid, ygrid, data = load_data(df, Rmid)
            N = 4096*4096
        except:
            xgrid, ygrid, data = load_snapshot_data(df, Rmid)
            N = 512*512
    spec = np.fft.rfftn(data)
    print(f"lambda = {path.Lambdaz[i]}, spec (min,max) = ({spec.min()},{spec.max()})")
    print(f"lambda = {path.Lambdaz[i]}, data (min,max) = ({data.min()},{data.max()})")
    im = ax.imshow(np.log10(np.abs(spec[:nk,:nm])**2/N + eps), vmin=-3,vmax=3, rasterized=True, origin='lower',interpolation='nearest')
    label = path.Lambdaz[i]
    if label is None:
        label = '\infty'
    ax.text(70,80,f"$\Lambda = {label}$", color='white')

# when cbar_mode is 'single', for ax in grid, ax.cax = grid.cbar_axes[0]

cbar = ax.cax.colorbar(im)
cbar = grid.cbar_axes[0].colorbar(im)
ax.cax.xaxis.set_label_position('top')
ax.cax.xaxis.set_ticks_position('top')
ax.cax.tick_params(labelsize=12)
cbar.set_label(r"$\log_{10} |\hat{v}|^2$")

plt.savefig("../figs/rei900_spectra.pdf", dpi=150)

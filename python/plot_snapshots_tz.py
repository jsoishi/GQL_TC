"""
Plot planes from joint analysis files.

Usage:
    plot_snapshots.py <files>... [--output=<dir> --period=<period>]

Options:
    --output=<dir>      Output directory [default: ./img_snapshots]
    --period=<period>   rotation period  [default: 1]
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import transforms
plt.ioff()
plt.style.use('prl')
import logging
logger = logging.getLogger(__name__)
from dedalus.extras import plot_tools


def main(filename, start, count, output, period):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    tasks = []
    with h5py.File(filename, mode='r') as file:
        keys = sorted(file['tasks'].keys())
        for task in keys:
            if "_r_slice" in task:
                tasks.append(task)
            
    scale = 2
    dpi = 300
    eta = 0.883
    Gamma = 29.9
    Rmid = eta/(1-eta) + 0.5
    Lx = Rmid*2*np.pi
    Lz = Gamma
    title_func = lambda sim_time: r'$t/\tau_\nu = {:.3f}$'.format(sim_time)
    savename_func = lambda write: 'urms_tz_{:06}.png'.format(write)
    # Layout
    nrows, ncols = 1,1
    image = plot_tools.Box(3*Lx/Lz, 3)
    pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
    margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)

    # Create multifigure
    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure
    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for index in range(start, start+count):
            logger.info(f"plotting index {index}")
            axes = mfig.add_axes(0, 0, [0, 0, 1, 1])
            # Call 3D plotting helper, slicing in time
            dset = np.sqrt(file['tasks']['u_r_slice'][index, :,:,0]**2 + file['tasks']['v_r_slice'][index, :,:,0]**2 + file['tasks']['w_r_slice'][index, :,:,0]**2)
            image_axes = [2, 1]
            data_slices = [index, slice(None), slice(None),0]
            xmesh, ymesh, vdata = plot_tools.get_plane(file['tasks/u_r_slice'], 2, 1, data_slices)
            xmesh *= Rmid

            # Setup axes
            # Bounds (left, bottom, width, height) relative-to-axes
            pbbox = transforms.Bbox.from_bounds(0.03, 0, 0.94, 0.94)
            cbbox = transforms.Bbox.from_bounds(0.03, 0.95, 0.94, 0.05)
            # Convert to relative-to-figure
            to_axes_bbox = transforms.BboxTransformTo(axes.get_position())
            pbbox = pbbox.transformed(to_axes_bbox)
            cbbox = cbbox.transformed(to_axes_bbox)
            # Create new axes and suppress base axes
            paxes = axes.figure.add_axes(pbbox)
            caxes = axes.figure.add_axes(cbbox)
            axes.axis('off')
            plot = paxes.pcolormesh(xmesh, ymesh, dset)
            paxes.axis(plot_tools.pad_limits(xmesh, ymesh))
            paxes.tick_params(length=0, width=0)

            # Colorbar
            cbar = plt.colorbar(plot, cax=caxes, orientation='horizontal',
                                ticks=ticker.MaxNLocator(nbins=5))
            cbar.outline.set_visible(False)
            caxes.xaxis.set_ticks_position('top')
            caxes.set_xlabel(r"$u'_{rms}$")
            caxes.xaxis.set_label_position('top')

            paxes.set_xlabel(r"$r_{mid} \theta$")
            paxes.set_ylabel(r"$z$")

            # Add time title
            title = title_func(file['scales/sim_time'][index]/period)
            title_height = 1 - 0.5 * mfig.margin.top / mfig.fig.y
            fig.suptitle(title, x=0.1, y=title_height, ha='left')
            # Save figure
            savename = savename_func(file['scales/write_number'][index])
            savepath = output.joinpath(savename)
            fig.savefig(str(savepath), dpi=dpi)
            fig.clear()
    plt.close(fig)


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = pathlib.Path(args['--output']).absolute()
    period = float(args['--period'])
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()
    post.visit_writes(args['<files>'], main, output=output_path, period=period)


"""
Plot planes from joint analysis files.

Usage:
    plot_snapshots.py <files>... [--output=<dir> --period=<period>]

Options:
    --output=<dir>  Output directory [default: ./img_snapshots]
    --period=<period>   rotation period  [default: 1]
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

from dedalus.extras import plot_tools


def main(filename, start, count, output, period):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    tasks = []
    with h5py.File(filename, mode='r') as file:
        keys = sorted(file['tasks'].keys())
        for task in keys:
            if "_t_slice" in task:
                tasks.append(task)
            
    scale = 4
    dpi = 100
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: 'write_rz_{:06}.png'.format(write)
    # Layout
    nrows, ncols = 1,3 
    image = plot_tools.Box(1, 3)
    pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
    margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)

    # Create multifigure
    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure
    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for index in range(start, start+count):
            for n, task in enumerate(tasks):
                # Build subfigure axes
                i, j = divmod(n, ncols)
                axes = mfig.add_axes(i, j, [0, 0, 1, 1])
                # Call 3D plotting helper, slicing in time
                dset = file['tasks'][task]
                image_axes = [3, 1]
                data_slices = [index, slice(None), 0, slice(None)]
                plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes, title=task, even_scale=True)
            # Add time title
            title = title_func(file['scales/sim_time'][index]/period)
            title_height = 1 - 0.5 * mfig.margin.top / mfig.fig.y
            fig.suptitle(title, x=0.48, y=title_height, ha='left')
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


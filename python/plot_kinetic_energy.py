import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import re
from timeseries import Timeseries

plt.style.use('prl')
filename = sys.argv[-1]
run_name = re.search("results/(.*)/scalar/.*",filename).group(1)
outfile = '../figs/KE_ens'+run_name+'.png'
ts = Timeseries(filename)
fig = ts.plot_en_ens()
for ax in fig.axes:
    ax.set_xlim(0,8)
    
fig.savefig(outfile)
print(f"final time: {ts.t[-1]/ts.period} inner rotation periods")

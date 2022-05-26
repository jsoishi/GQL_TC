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
fig.savefig(outfile)

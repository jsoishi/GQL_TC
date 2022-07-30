import re
import h5py
import numpy as np
from pathlib import Path

class RunSeries:
    """Tool to manage a series of runs.

    """
    def __init__(self, run_list, name=None, root_dir="results/"):
        self.run_list = run_list
        self.name = name
        self.root_dir = Path(root_dir)
        self.sim_time = []
        self.data = {}

        self.start_time = []
        self.eta = []
        self.mu = []
        self.Re = []
        self.Re_o = []
        self.Lz = []
        self.Lambdaz = []
        self.Lambdat = []

        for r in self.run_list:
            eta, Re, Lz, mu, R1, R2, Omega1, period, Re_o, Lambdaz, Lambdat = self.parse_run_name(r)
            self.eta.append(eta)
            self.Lz.append(Lz)
            self.mu.append(mu)
            self.Re.append(Re)
            self.Re_o.append(Re_o)
            self.Lambdaz.append(Lambdaz)
            self.Lambdat.append(Lambdat)

    def parse_run_name(self, fn):
        
        eta = float(re.search("eta_([\d.e+-]+)",fn).group(1))
        Re = float(re.search("re_([\d.e+-]+)",fn).group(1))
        Lz = float(re.search("Gamma_([\d.e+-]+)",fn).group(1))
        mu = float(re.search("mu_([\d.e+-]+)",fn).group(1))
        Lambdaz = re.search("Lambdaz\_(\d+)",fn)
        if Lambdaz:
            Lambdaz = int(Lambdaz.group(1))
        else:
            Lambdaz = None
        Lambdat = re.search("Lambdat\_(\d+)",fn)
        if Lambdat:
            Lambdat = int(Lambdat.group(1))
        else:
            Lambdat = None

        # derived parameters
        R1 = eta/(1. - eta)
        R2 = 1./(1-eta)
        Omega1 = 1/R1
        period = 2*np.pi/Omega1
        Re_o = Re*mu/eta

        return eta, Re, Lz, mu, R1, R2, Omega1, period, Re_o, Lambdaz, Lambdat
        
    def load(self, filetype, field, concatenate=True):
        filetype = Path(filetype)
        self.data[field] = []
        for r in self.run_list:
            # make sure there's only one 
            filename = self.root_dir / Path(r) / filetype / f"{filetype}.h5"
            with h5py.File(filename, "r") as df:
                self.sim_time.append(df['scales/sim_time'][:])
                self.data[field].append(df[f'tasks/{field}'][:])
            self.start_time.append(self.sim_time[-1][0])
        if concatenate:
            self.sim_time = np.concatenate(self.sim_time)
            self.data[field] = np.concatenate(self.data[field])

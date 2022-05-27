import h5py
import numpy as np
import re
import matplotlib.pyplot as plt
class Timeseries:
    def __init__(self, filename, ref_growth=None):
        self.fn = filename
        self.ref_growth = ref_growth
        
        self.eta = float(re.search("eta_([\d.e+-]+)",self.fn).group(1))
        self.R1 = self.eta/(1. - self.eta)
        self.R2 = 1./(1-self.eta)
        self.Omega1 = 1/self.R1
        self.period = 2*np.pi/self.Omega1
        self.Re = float(re.search("re_([\d.e+-]+)",self.fn).group(1))
        self.Lz = float(re.search("Gamma_([\d.e+-]+)",self.fn).group(1))
        with h5py.File(self.fn, "r") as ts:
            self.t = ts['scales/sim_time'][:]
            self.w_rms = ts['tasks/w_rms'][:,0,0,0]
            self.KE = ts['tasks/KE'][:,0,0,0]/self.Lz
            self.KEp = ts['tasks/pertubation_KE'][:,0,0,0]/self.Lz
            self.u_p = ts['tasks/u_probe'][:,0,0,0]
            self.w_p = ts['tasks/w_probe'][:,0,0,0]
            self.enstrophy = ts['tasks/enstrophy'][:,0,0,0]/self.Lz

        
    def analyze_growth(self, window=(2,14)):
        self.t_window = (self.t/self.period > window[0]) & (self.t/self.period < window[1])
        self.gamma_w, self.log_w0 = np.polyfit(self.t[self.t_window], np.log(self.w_rms[self.t_window]),1)

    @property
    def rel_error(self):
        if self.ref_growth is None:
            raise ValueError("Unspecified reference value for relative error calculation.")

        rel_error = (self.gamma_w - self.ref_growth)/self.ref_growth
        return rel_error

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.semilogy(self.t/self.period,self.w_rms)
        ax.semilogy(self.t/self.period, np.exp(self.log_w0)*np.exp(self.gamma_w*self.t), '-.', label='$\gamma_w = %f$' % self.gamma_w)
        ax.set_ylim(np.min(self.w_rms), np.max(3*self.w_rms))
        ax.set_xlabel(r"$t/\tau_{1}$")
        ax.set_ylabel(r"$w_{rms}$")
        ax.legend()
        plt.tight_layout()

        return ax
        
    def plot_en_ens(self, fig=None):

        if fig is None:
            fig = plt.figure()
        KE_ax = fig.add_subplot(211)
        en_ax = fig.add_subplot(212)
        u_to_visc = lambda x: x/self.Re*self.period
        visc_to_u = lambda x: self.Re*x/self.period
        ax2 = KE_ax.secondary_xaxis("top",functions=(u_to_visc, visc_to_u))
        ax2.set_xlabel(r"$t/\tau_{\nu}$")
        #KE_ax.plot(self.t/self.period, self.KE, label='total')
        KE_ax.plot(self.t/self.period, self.KEp, label='perturbation')
        KE_ax.set_xlabel(r"$t/\tau_{1}$")
        KE_ax.set_ylabel(r"$KE/L_z$")

        en_ax.plot(self.t/self.period, self.enstrophy)
        en_ax.set_xlabel(r"$t/\tau_{1}$")
        en_ax.set_ylabel(r"$\mathcal{E}/L_z$")
        plt.tight_layout()

        return fig

        

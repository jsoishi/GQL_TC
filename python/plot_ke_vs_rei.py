from series import RunSeries
import matplotlib.pyplot as plt
import numpy as np
from paths import path2, path1
mes1 = np.loadtxt("meseguer_fig_4_path1.dat", delimiter=',')
mes2 = np.loadtxt("meseguer_fig_4_path2.dat", delimiter=',')
plt.style.use('prl')
field = "pertubation_KE"
path2 = RunSeries(path2)
path2.load("scalar", field, concatenate=False)

path1 = RunSeries(path1)
path1.load("scalar", field, concatenate=False)
KE_rms_2 = []
eta = path1.eta[0]
R2 = 1/(1-eta)
R1 = eta/(1-eta)
Vol = np.pi*path1.Lz[0]*(R2**2 - R1**2)
for i, Re in enumerate(path2.Re):
    KE_rms_2.append(np.sqrt((path2.data[field][i]/Vol).mean()))
KE_rms_1 = []
for i, Re in enumerate(path1.Re):
    KE_rms_1.append(np.sqrt((path1.data[field][i]/Vol).mean()))
    
plt.scatter(path1.Re, KE_rms_1, label=r'$Re_o=-3398$')
plt.scatter(path2.Re, KE_rms_2, label=r'$Re_o=-1359$')
plt.scatter(mes1[:,0], mes1[:,1]/mes1[:,0],label='Path 1 Meseguer et al')
plt.scatter(mes2[:,0], mes2[:,1]/mes2[:,0],label='Path 2 Meseguer et al')
plt.legend()
plt.xlabel(r"$\mathrm{Re}_i$")
plt.ylabel(r"$\left<KE_{pert}\right>$")
plt.tight_layout()
plt.savefig(f'../figs/ke_vs_rei_our_units.png', dpi=300)

"""
This script plots the angular momentum of the flow at any given time
and compares the results to the Marcus 2 paper. Could easily be edited
to add more radial variables to plot.

python3 plot_radial_data.py folder_prefix

"""
import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as de

import os
import glob

# Parses filename passed to script
parser = argparse.ArgumentParser(description='Unique folder prefix to identify all runs to plot')
parser.add_argument('folder_prefix', metavar='Rc', type=str, help='.h5 file to plot radial data')
args = parser.parse_args()
prefix = vars(args)['folder_prefix']

prefix_ast = str(prefix) + "*"

eta = 0.875
mu = 0.

# define couette flow
def couette(r):
    A = (1/eta - 1.)*(mu-eta**2)/(1-eta**2)
    B = eta*(1-mu)/((1-eta)*(1-eta**2))
    v0 = A*r + B/r
    return v0

# GQL_arr = []
# DNS_arr = []
# lambda_list = []

folders = sorted(glob.glob(prefix_ast), key=lambda folder: int(folder.split("_")[-1]))
if 'GQL' not in folders[-1]:
    folders.insert(0,folders[-1])
    del folders[-1]
else:
    folders.insert(0,folders[-2])
    del folders[-2]

for folder in folders:
    slices = folder + "/slices/slices_s1.h5"
    datafile = h5py.File(slices,'r')
    r_vs_v_plane_avg = datafile['tasks/v_tot'][-100:,0,0,:].mean(axis=0)
    r = datafile['scales/r/1.0'][:]
    if "GQL" not in folder:
        plt.plot(r,r_vs_v_plane_avg,label="DNS", linewidth=3)
        # DNS_arr.append(r_vs_v_plane_avg)
    else:
        plt.plot(r,r_vs_v_plane_avg,label="GQL, $\Lambda = " + folder.split("_")[-1] + "$")
        # GQL_arr.append([folder.split("_")[-1], r_vs_v_plane_avg])
        # lambda_list.append(int(folder.split("_")[-1]))

plt.legend()

plt.xlabel("$r$")
plt.ylabel("Avg. V velocity")
plt.tight_layout()
plt.style.use('prl')
plt.show()
plot_file_name = 'radial_profile_multiple_' + str(prefix)[0:-1] + '.png'
plt.savefig(plot_file_name, dpi=300)

# GQL_arr = np.array(GQL_arr)
# DNS_arr = np.array(DNS_arr)

# nr = 32
# R1 = 7
# R2 = 8

# r = de.Chebyshev('r', nr, interval=(R1,R2))
# d = de.Domain([r], grid_dtype=np.float)
# v_r_DNS = d.new_field()
# v_r_GQL = d.new_field()
# r = d.grid(0)

# lambda_vals = np.array(lambda_list)
# l2_vals = np.zeros(len(lambda_vals))

# for i in range(len(GQL_arr)):
#     v_r_DNS['g'] = DNS_arr[0]
#     v_r_GQL['g'] = GQL_arr[i][1]
#     derived_field = ((v_r_GQL - v_r_DNS)**2).evaluate() # must call evaluate()
#     integ_derived_field = derived_field.integrate()
#     l2_vals[i] = integ_derived_field['g'][0]

# plt.plot(lambda_vals, l2_vals, label = 'L2 norm ' + str(prefix)[0:-1], marker = 'o')

# plt.legend()

# plt.title("L2 norm vs. $\Lambda$")
# plt.xlabel("$\Lambda$")
# plt.ylabel("L2 norm")
# plt.tight_layout()
# plt.style.use('prl')
# plt.show()

# plot_file_name = 'L2norm_' + str(prefix) + '.png'
# plt.savefig(plot_file_name, dpi=300)

# plt.semilogy(lambda_vals, l2_vals, label = 'L2 norm', marker = 'o')

# plt.legend()

# plt.title("L2 norm vs. $\Lambda$")
# plt.xlabel("$\Lambda$")
# plt.ylabel("L2 norm")
# plt.tight_layout()
# plt.style.use('prl')
# plt.show()

# plot_file_name = 'L2norm_semilog_' + str(prefix)[0:-1] + '.png'
# plt.savefig(plot_file_name, dpi=300)
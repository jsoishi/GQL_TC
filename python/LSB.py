import matplotlib.pyplot as plt
import numpy as np

plt.style.use('prl')
andereck = np.loadtxt("data/andereck_lsb.csv",delimiter=',')
ours     = np.loadtxt("data/our_LSB.csv",delimiter=',')
our_path=  np.loadtxt('data/our_path_1.csv',delimiter=',')
our_path_2=np.loadtxt('data/our_path_2.csv',delimiter=',')

# conversion factors for Coles (1965)
# ri_c = 5.
# ro_c = 11.44/2
# d_c = ro_c - ri_c
# ifactor = d_c/ri_c
# ofactor = d_c/ro_c

# coles_H1z = np.loadtxt("data/coles_1965_zoom.csv",delimiter=',')
# cz_x = ofactor*(coles_H1z[:,2] + coles_H1z[:,0])/2
# cz_y= ifactor*(coles_H1z[:,3] + coles_H1z[:,1])/2
# cz_xerr = ofactor*(coles_H1z[:,2] - coles_H1z[:,0])/2
# cz_yerr = ifactor*(coles_H1z[:,3] - coles_H1z[:,1])/2

# andereck = np.loadtxt("data/andereck_H2.csv",delimiter=',')
# andereck_x    = (andereck[:,2] + andereck[:,0])/2
# andereck_y    = (andereck[:,3] + andereck[:,1])/2
# andereck_xerr = (andereck[:,2] - andereck[:,0])/2
# andereck_yerr = (andereck[:,3] - andereck[:,1])/2

# coles_H1 = np.loadtxt("data/coles_1965_H1.csv",delimiter=',')
# c_xerr = ofactor*(coles_H1[:,2] - coles_H1[:,0])
# c_yerr = ifactor*(coles_H1[:,3] - coles_H1[:,1])

# meseguer_H1 = np.loadtxt("data/meseguer_H1.csv",delimiter=',')
# meseguer_H2 = np.loadtxt("data/meseguer_H2.csv",delimiter=',')
# deguchi_H1 = np.loadtxt("data/deguchi_etal_coles.csv", delimiter=',')
# d_xerr = deguchi_H1[:,2] - deguchi_H1[:,0]
# d_yerr = deguchi_H1[:,3] - deguchi_H1[:,1]
# ad_H2 = np.loadtxt("data/andereck_H2.csv",delimiter=',')

#plt.plot(andereck[:,0], andereck[:,1], 'x-',label='Andereck et al (1986)')
#plt.plot(meseguer[:,0], meseguer[:,1], '*-',label='Meseguer et al (2009)')
#plt.plot(ours[:,0], ours[:,1], '^-',label='ours')
# plt.errorbar(coles_H1[:,0]*ofactor, coles_H1[:,1]*ifactor, xerr=c_xerr,yerr=c_yerr, label='Coles (1965) 2a')
# plt.errorbar(cz_x, cz_y, xerr=cz_xerr,yerr=cz_yerr, label='Coles (1965) 2b')
# plt.errorbar(andereck_x, andereck_y, xerr=andereck_xerr,yerr=andereck_yerr, label='Andereck et al (1986)')

#plt.plot(coles_H2[:,0]*ofactor, coles_H2[:,1]*ifactor, 'o', label='Coles (1965)')
#plt.plot(our_path[:,0],our_path[:,1],'k+')
#plt.plot(our_path_2[:,0],our_path_2[:,1],'k+')

# plt.plot(meseguer_H1[:,0], meseguer_H1[:,1],'^', label='Meseguer et al (2009) H1')
# plt.errorbar(deguchi_H1[:,0], deguchi_H1[:,1],xerr=d_xerr, yerr=d_yerr, label='Deguchi et al (2014)')
# plt.plot(meseguer_H2[:,0], meseguer_H2[:,1],'^', label='Meseguer et al (2009) H2')
# plt.legend(loc='lower right',fontsize=14)
# plt.xlabel(r"$\mathrm{Re}_o$")
# plt.ylabel(r"$\mathrm{Re}_i$")



Rei_g = [861.53846154, 810.25641026, 758.97435897, 707.69230769, 656.41025641, 605.12820513, 553.84615385, 502.56410256, 451.28205128, 400.]
Reo_g = np.linspace(-4000,-1000, 10)
plt.figure()
plt.plot(ours[:,0], ours[:,1])
plt.plot(our_path[:3,0],our_path[:3,1],'k+',alpha=0.3, label='DNS only')
plt.plot(our_path[:3:2,0],our_path[:3:2,1],'ko', label='DNS & GQL')
plt.plot(our_path[:3,0],our_path[:3,1],'k',alpha=0.3,zorder=1,lw=1)
plt.annotate("", xy=(our_path[0,0],750), xytext=(our_path[0,0], 760),
             arrowprops=dict(arrowstyle="-|>",facecolor='k',edgecolor='None',alpha=0.4),zorder=3)

plt.plot(our_path_2[:-2,0],our_path_2[:-2,1],'k+', alpha=0.3)
plt.plot(our_path_2[:-2,0],our_path_2[:-2,1],'k',alpha=0.3,lw=1)
#plt.plot((our_path_2[0,0],our_path_2[0,0]),(our_path_2[0,1], our_path_2[6,1]),'k+')
plt.plot((our_path_2[0,0],our_path_2[0,0]),(our_path_2[0,1], our_path_2[6,1]),'ko')
plt.annotate("", xy=(our_path_2[0,0],570), xytext=(our_path_2[0,0], 580),
             arrowprops=dict(arrowstyle="-|>",facecolor='k',edgecolor='None',alpha=0.4))

#plt.legend(loc='lower right',fontsize=14)
plt.xlim(-3500,-1000)
plt.xlabel(r"$\mathrm{Re}_o$")
plt.ylabel(r"$\mathrm{Re}_i$")
plt.legend()
plt.tight_layout()
plt.savefig("../figs/reo_rei_lsb.pdf")

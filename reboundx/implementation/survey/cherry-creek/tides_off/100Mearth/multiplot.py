import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# load REBOUND data
data = np.loadtxt('output/m.txt')  # return (N, 2) array
ts = data[:, 0]                    # return only 1st col
mass = data[:, 1]                  # return only 2nd col
data = np.loadtxt('output/r.txt')
radius = data[:, 1]                # data in AU
init_as = np.arange(0.4, 1.51, 0.1) # in AU
a = np.zeros([ts.size, init_as.size])
for i,init_a in enumerate(init_as):
    fname = 'output/{:.1f}au.txt'.format(init_a)
    data = np.loadtxt(fname)
    a[:, i] = data[:, 1]           # data in AU

f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
# fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8, 12), dpi=300)
fig.subplots_adjust(hspace=0)

ax1.set_ylabel("$M(t)$ / $M_{\odot}$", fontsize='large')
ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator())
ax1.plot(ts,mass, color='black', lw=2.5)
ax1.grid()

ax2.set_xlabel('Time / yr', fontsize='large')
ax2.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator())
ax2.set_ylabel('Distance / AU', fontsize='large')
ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator())
ax2.plot(ts,radius, color='black', lw=2.5, label='$R(t)$')
for i,init_a in enumerate(init_as):
    # label = '$a_{:.1f}(t)$'%(init_a)
    ax2.plot(ts,a[:, i])
ax2.legend(fontsize='large', loc='best')
ax2.grid()

plt.gca().xaxis.set_major_formatter(mticker.FuncFormatter(g))
plt.show()
# plt.savefig('survey.png')

import numpy as np

from astropy.io import ascii

import matplotlib.pyplot as plt
#plt.style.use(style.notebook)

import plotsettings as plotsettings
import matplotlib as matplotlib
#matplotlib.use('pdf')

publishable = plotsettings.Set('Cell')
publishable.set_figsize(1.1, 1.0, aspect_ratio=0.95)
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
matplotlib.rc('text', usetex=False)

files = ["lognew.txt", "lognew2.txt"]

for file in files:
    print(file)
    data = ascii.read(file, data_start=1)
    print data
    kin = data['col2']
    pot   = data['col4']
    tot = -pot+kin

    #print len(kin)
    #plt.plot(np.arange(0, len(pot)), kin)
    #plt.plot(np.arange(0, len(pot)), pot)

    plt.plot(np.arange(0, len(pot))*0.1, np.abs(tot/tot[0]-1), lw=0.8)
    #plt.plot(P, time1[0]*0.08*np.power(5.5**P, 2.5), 'b--', alpha=0.8, lw=2.)

plt.xlabel(r"$time $")
plt.ylabel(r"$energy error$")

#plt.xlim(1, 10)
#plt.ylim(-3, 3)
#plt.xscale("log")
#plt.yscale("log")
#plt.legend(loc='upper left', ncol=1, fontsize='x-small')
plt.yscale("log")
plt.grid(True)
plt.tight_layout()
plt.savefig("energy_error.pdf")
plt.close()

